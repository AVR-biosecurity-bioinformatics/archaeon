process GFFREAD_EXTRACT {

    tag 'representative receptor sequences'

    conda "${moduleDir}/environment.yml"

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