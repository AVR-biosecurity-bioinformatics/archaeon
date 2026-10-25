process FILTER_CANDIDATES {

    tag 'filter candidate models'

    conda "${moduleDir}/environment.yml"

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