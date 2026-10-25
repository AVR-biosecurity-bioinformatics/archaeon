process SELECT_REPRESENTATIVE_LOCI {

    tag 'representative receptor loci'
    conda "${moduleDir}/environment.yml"

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