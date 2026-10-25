params {

    // Mandatory inputs
    genome: Path
    references: Path
    reference_metadata: Path
    slurm_account: String                      // SLURM account for submitted jobs.

    // Optional inputs
    min_coverage: Float = 0.5f
    min_identity: Float = 0.3f
    max_intron: Integer = 200000
    max_hits: Integer = 50
    secondary_ratio: Float = 0.5f
    translation_table: Integer = 1

    locus_cluster_distance: Integer = 0

}

/*
    Define Nextflow processes
    Each process wraps a tool or a discrete processing step
*/

process MINIPROT_INDEX {
    tag 'genome'
    conda 'bioconda::miniprot=0.18'
    cpus 4
    memory '8 GB'

    input:
    path genome, stageAs: 'genome.fa'

    output:
    path 'genome.mpi', emit: index

    script:
    """
    miniprot -t ${task.cpus} -T ${params.translation_table} -d genome.mpi genome.fa
    """
}

process MINIPROT_ALIGN {

    tag 'reference panel'

    conda 'bioconda::miniprot=0.18 conda-forge::gawk'

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

// Retain protein-to-genome alignments passing the minimum query-coverage
// and alignment-identity thresholds.

process FILTER_CANDIDATES {

    tag 'filter candidate models'

    conda 'conda-forge::gawk'

    cpus 1
    memory '1 GB'

    input:
    path models, stageAs: 'miniprot.models.tsv'
    path gff, stageAs: 'miniprot.gff3'

    output:
    path 'candidates.tsv', emit: summary
    path 'candidates.gff3', emit: gff

    script:
    """
    # Retain models passing minimum coverage and identity thresholds.
    gawk -F '\\t' \
        -v min_coverage=${params.min_coverage} \
        -v min_identity=${params.min_identity} \
        'NR == 1 || (\$8 >= min_coverage && \$9 >= min_identity)' \
        miniprot.models.tsv \
        > candidates.tsv

    # Retain candidate mRNA models and their associated child features.
    gawk -F '\\t' '
    FNR == NR {
        if (FNR > 1)
            keep[\$1] = 1

        next
    }

    /^##gff-version/ {
        if (!header++)
            print "##gff-version 3"

        next
    }

    /^#/ {
        next
    }

    {
        delete id
        delete parent

        match(\$9, /(^|;)ID=([^;]+)/, id)
        match(\$9, /(^|;)Parent=([^;]+)/, parent)

        if ((\$3 == "mRNA" && id[2] in keep) || (\$3 != "mRNA" && parent[2] in keep))
            print
    }
    ' candidates.tsv miniprot.gff3 \
        > candidates.gff3
    """
}

/*
    Group overlapping, same-strand candidate models into genomic loci
    and select the highest-supported model from each locus.

    Representative models are ranked by:
      1. alignment score, descending
      2. identity, descending
      3. query coverage, descending
      4. model ID, ascending
*/
process SELECT_REPRESENTATIVE_LOCI {

    tag 'representative receptor loci'

    conda 'bioconda::bedtools=2.31.1 conda-forge::gawk'

    cpus 1
    memory '1 GB'

    input:
    path candidates, stageAs: 'candidates.tsv'

    output:
    path 'representative_models.tsv', emit: models

    script:
    """
    # Convert candidate models to BED, cluster overlapping same-strand
    # intervals, rank models within each cluster, and retain the best model.
    {
        # Write a tab-delimited header. The escapes are part of the printf
        # format string, so they are converted into real tab characters.
        printf 'locus_id\\tmodel_id\\tquery_id\\tcontig\\tstart\\tend\\tstrand\\tquery_coverage\\tidentity\\talignment_score\\n'

        awk -F '\\t' '
        BEGIN {
            OFS = "\\t"
        }

        NR > 1 {
            # Convert one-based inclusive coordinates to BED coordinates.
            print \
                \$4, \
                \$5 - 1, \
                \$6, \
                \$1, \
                \$10, \
                \$7, \
                \$2, \
                \$8, \
                \$9
        }
        ' candidates.tsv \
        | LC_ALL=C sort \
            -k1,1 \
            -k2,2n \
            -k3,3n \
            -k6,6 \
        | bedtools cluster \
            -s \
            -d ${params.locus_cluster_distance} \
            -i - \
        | LC_ALL=C sort \
            -k10,10n \
            -k5,5nr \
            -k9,9nr \
            -k8,8nr \
            -k4,4 \
        | awk -F '\\t' '
        BEGIN {
            OFS = "\\t"
        }

        # After sorting, the first model in each cluster is the
        # provisional representative for that genomic locus.
        !seen[\$10]++ {
            print \
                sprintf("LOCUS_%06d", \$10), \
                \$4, \
                \$7, \
                \$1, \
                \$2 + 1, \
                \$3, \
                \$6, \
                \$8, \
                \$9, \
                \$5
        }
        '
    } > representative_models.tsv
    """
}


/*
    Extract spliced CDS and translated protein sequences for the provisional
    representative model selected from each candidate genomic locus.

    Representative model IDs are used to subset the candidate GFF3 before
    sequence extraction.

    -x writes spliced CDS nucleotide sequences.
    -y writes translated protein sequences.
*/
process GFFREAD_EXTRACT {

    tag 'representative receptor sequences'

    conda 'bioconda::gffread=0.12.9 conda-forge::gawk'

    cpus 1
    memory '4 GB'

    input:
    path genome, stageAs: 'genome.fa'
    path representatives, stageAs: 'representative_models.tsv'
    path candidate_gff, stageAs: 'candidates.gff3'

    output:
    path 'representative_models.cds.fna', emit: cds
    path 'representative_models.proteins.faa', emit: proteins

    script:
    """
    # Subset the candidate GFF3 to the selected representative mRNAs
    # and their associated CDS and stop-codon features.
    gawk -F '\\t' '
    BEGIN {
        OFS = "\\t"
    }

    # Load representative model IDs from the first input.
    FNR == NR {
        if (FNR > 1)
            keep[\$2] = 1

        next
    }

    # Write one clean GFF3 header and discard other comments.
    /^#/ {
        if (!header++)
            print "##gff-version 3"

        next
    }

    {
        delete attr

        n = split(\$9, fields, ";")

        for (i = 1; i <= n; i++) {
            split(fields[i], kv, "=")
            attr[kv[1]] = kv[2]
        }

        # Retain selected mRNA records.
        if (\$3 == "mRNA") {
            if (attr["ID"] in keep)
                print

            next
        }

        # Retain child features belonging to selected mRNAs.
        n = split(attr["Parent"], parents, ",")

        for (i = 1; i <= n; i++) {
            if (parents[i] in keep) {
                print
                break
            }
        }
    }
    ' representative_models.tsv candidates.gff3 \
    > representative_models.gff3

    # Extract the spliced CDS and translated protein for each
    # representative gene model.
    gffread \
        representative_models.gff3 \
        -g genome.fa \
        -x representative_models.cds.fna \
        -y representative_models.proteins.faa

    test -s representative_models.cds.fna
    test -s representative_models.proteins.faa
    """
}

/*
    Assign each representative locus to the nAChR subunit associated with
    its highest-supported reference protein.

    These are provisional similarity-based assignments, not confirmed
    orthology calls.

    This should be updated using phylogenetic classificaiton later on
*/
process CLASSIFY_REPRESENTATIVE_LOCI {

    tag 'classify receptor loci'

    conda 'conda-forge::gawk'

    cpus 1
    memory '1 GB'

    input:
    path representatives, stageAs: 'representative_models.tsv'
    path metadata, stageAs: 'reference_metadata.tsv'

    output:
    path 'classified_receptors.tsv', emit: classifications

    script:
    """
    gawk -F '\\t' '
    BEGIN {
        OFS = "\\t"
    }

    # Load reference accession metadata.
    FNR == NR {
        if (FNR > 1) {
            species[\$1] = \$2
            subunit[\$1] = \$3
        }

        next
    }

    # Add subunit classification and confidence fields.
    FNR == 1 {
        print \
            \$0, \
            "reference_species", \
            "subunit", \
            "assignment_confidence"

        next
    }

    {
        query = \$3
        coverage = \$8
        identity = \$9

        if (!(query in subunit)) {
            print \
                "ERROR: no metadata found for reference " query \
                > "/dev/stderr"

            exit 1
        }

        if (coverage >= 0.90 && identity >= 0.80)
            confidence = "high"
        else if (coverage >= 0.80 && identity >= 0.60)
            confidence = "moderate"
        else
            confidence = "low"

        print \
            \$0, \
            species[query], \
            subunit[query], \
            confidence
    }
    ' reference_metadata.tsv representative_models.tsv \
    > classified_receptors.tsv
    """
}
/*
    Nextflow workflow 
    This is where processes are connected and pass their outputs to downstream steps
*/

workflow {

    main:

    // Create value channels for the target genome and reference proteins.
    genome = channel.value(params.genome)
    references = channel.value(params.references)
    metadata = channel.value(params.reference_metadata)

    // Build a searchable Miniprot index of the target genome.
    MINIPROT_INDEX(
        genome
    )

    // Align reference receptor proteins directly to the genome and produce:
    //   1. Miniprot GFF3 gene models
    //   2. PAF alignment records
    //   3. A compact model summary table
    MINIPROT_ALIGN(
        MINIPROT_INDEX.out.index,
        references
    )

    // Retain alignments passing the minimum query-coverage and identity
    // thresholds, and subset the Miniprot GFF3 to those candidate models.
    FILTER_CANDIDATES(
        MINIPROT_ALIGN.out.models,
        MINIPROT_ALIGN.out.gff
    )

    // Cluster overlapping, same-strand models and retain the
    // highest-supported model from each putative genomic locus.
    SELECT_REPRESENTATIVE_LOCI(
        FILTER_CANDIDATES.out.summary
    )

    // Assign each representative locus using the reference protein that
    // produced the highest-supported model at that genomic locus.
    CLASSIFY_REPRESENTATIVE_LOCI(
        SELECT_REPRESENTATIVE_LOCI.out.models,
        metadata
    )
    
    // Extract CDS and protein sequences for each representative model.
    // This can run in parallel with classification because both processes
    // depend only on the selected representative models.
    GFFREAD_EXTRACT(
        genome,
        SELECT_REPRESENTATIVE_LOCI.out.models,
        FILTER_CANDIDATES.out.gff
    )

    publish:

    // Raw Miniprot alignment outputs.
    miniprot_gff = MINIPROT_ALIGN.out.gff
    miniprot_paf = MINIPROT_ALIGN.out.paf
    miniprot_models = MINIPROT_ALIGN.out.models

    // Models passing the coverage and identity thresholds.
    candidates_tsv = FILTER_CANDIDATES.out.summary
    candidates_gff = FILTER_CANDIDATES.out.gff

    // Provisional representative model selected from each genomic locus.
    representative_models = SELECT_REPRESENTATIVE_LOCI.out.models

    // CDS and protein sequences extracted from representative models.
    representative_cds = GFFREAD_EXTRACT.out.cds
    representative_proteins = GFFREAD_EXTRACT.out.proteins

    // Classified receptors
    classified_receptors = CLASSIFY_REPRESENTATIVE_LOCI.out.classifications

}

/*
    Define the files published from the entry workflow.
*/
output {

    miniprot_gff {
        path 'alignments'
    }

    miniprot_paf {
        path 'alignments'
    }

    miniprot_models {
        path 'alignments'
    }

    candidates_tsv {
        path 'candidates'
    }

    candidates_gff {
        path 'candidates'
    }

    representative_models {
        path 'representatives'
    }

    representative_cds {
        path 'representatives'
    }

    representative_proteins {
        path 'representatives'
    }

    classified_receptors {
        path 'classification'
    }
}