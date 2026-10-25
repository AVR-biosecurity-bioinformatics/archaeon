process CLASSIFY_REPRESENTATIVE_LOCI {

    tag 'classify receptor loci'
    conda "${moduleDir}/environment.yml"

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