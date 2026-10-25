process MINIPROT_ALIGN {

    tag 'reference panel'

    conda "${moduleDir}/environment.yml"

    cpus 4
    memory '8 GB'

    input:
    path index, stageAs: 'genome.mpi'
    path references, stageAs: 'references.faa'

    output:
    path 'miniprot.gff3', emit: gff
    path 'miniprot.paf', emit: paf
    path 'miniprot.models.tsv', emit: models

    script:
    """
    # Align reference proteins to the target genome. --gff emits GFF3
    # features together with embedded ##PAF alignment records.
    miniprot \
        -t ${task.cpus} \
        -j 2 \
        -G ${params.max_intron} \
        -N ${params.max_hits} \
        --outn ${params.max_hits} \
        -p ${params.secondary_ratio} \
        --outs ${params.secondary_ratio} \
        --outc ${params.min_coverage} \
        --gff \
        genome.mpi \
        references.faa \
        > miniprot.gff3

    # Extract the embedded PAF records.
    awk '
        /^##PAF/ {
            sub(/^##PAF[[:space:]]*/, "")
            print
        }
    ' miniprot.gff3 > miniprot.paf

   # Create a compact model table for filtering and locus clustering.
    printf '%s\n' \
    'model_id\tquery_id\tquery_rank\tcontig\tstart\tend\tstrand\tquery_coverage\tidentity\talignment_score' \
    > miniprot.models.tsv

    # Join PAF-derived query coverage to GFF3 models using query ID and rank.
    awk -F '\t' '
    BEGIN {
        OFS = "\t"
    }

    # Load coverage from the embedded PAF records.
    FNR == NR {
        rank = ++query_rank[\$1]
        coverage[\$1 SUBSEP rank] = (\$4 - \$3) / \$2
        next
    }

    # Summarise each GFF3 mRNA model.
    \$3 == "mRNA" {
        delete attr

        n = split(\$9, fields, ";")

        for (i = 1; i <= n; i++) {
            split(fields[i], kv, "=")
            attr[kv[1]] = substr(fields[i], length(kv[1]) + 2)
        }

        # Target contains: query_id query_start query_end.
        split(attr["Target"], target, " ")

        query = target[1]
        rank = attr["Rank"]
        key = query SUBSEP rank

        if (!(key in coverage)) {
            print \
                "ERROR: no PAF alignment for query " query \
                ", rank " rank \
                ", model " attr["ID"] \
                > "/dev/stderr"

            exit 1
        }

        print \
            attr["ID"], \
            query, \
            rank, \
            \$1, \
            \$4, \
            \$5, \
            \$7, \
            sprintf("%.6f", coverage[key]), \
            attr["Identity"], \
            \$6
    }
    ' miniprot.paf miniprot.gff3 \
    >> miniprot.models.tsv

    """
}