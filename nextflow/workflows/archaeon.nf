/*
    Main Archaeon workflow.

    Identifies candidate receptor loci, selects representative gene
    models, assigns provisional receptor classifications, and extracts
    representative CDS and protein sequences.
*/

include { IDENTIFY_LOCI } from '../subworkflows/identify_loci'


workflow ARCHAEON {

    take:

    genome
    references
    metadata

    main:

    IDENTIFY_LOCI(
        genome,
        references,
        metadata
    )

    emit:

    // Raw Miniprot alignment outputs.
    miniprot_gff = IDENTIFY_LOCI.out.miniprot_gff
    miniprot_paf = IDENTIFY_LOCI.out.miniprot_paf
    miniprot_models = IDENTIFY_LOCI.out.miniprot_models

    // Models passing the coverage and identity thresholds.
    candidates_tsv = IDENTIFY_LOCI.out.candidates_tsv
    candidates_gff = IDENTIFY_LOCI.out.candidates_gff

    // Representative model selected from each genomic locus.
    representative_models = IDENTIFY_LOCI.out.representative_models

    // Extracted representative sequences.
    representative_cds = IDENTIFY_LOCI.out.representative_cds
    representative_proteins = IDENTIFY_LOCI.out.representative_proteins

    // Best-reference receptor classifications.
    classified_receptors = IDENTIFY_LOCI.out.classified_receptors
}