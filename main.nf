params {

    // Mandatory inputs
    genome: Path
    references: Path
    reference_metadata: Path
    slurm_account: String

    // Optional inputs
    min_coverage: Float = 0.5f
    min_identity: Float = 0.3f
    max_intron: Integer = 200000
    max_hits: Integer = 50
    secondary_ratio: Float = 0.5f
    translation_table: Integer = 1
    locus_cluster_distance: Integer = 0
}


include { ARCHAEON } from './nextflow/workflows/archaeon'


/*
    Archaeon entry workflow.
*/
workflow {

    main:

    // Create value channels for the target genome, reference proteins,
    // and reference metadata.
    genome = channel.value(params.genome)
    references = channel.value(params.references)
    metadata = channel.value(params.reference_metadata)

    // Run the complete Archaeon workflow.
    ARCHAEON(
        genome,
        references,
        metadata
    )

    publish:

    // Raw Miniprot alignment outputs.
    miniprot_gff = ARCHAEON.out.miniprot_gff
    miniprot_paf = ARCHAEON.out.miniprot_paf
    miniprot_models = ARCHAEON.out.miniprot_models

    // Models passing the coverage and identity thresholds.
    candidates_tsv = ARCHAEON.out.candidates_tsv
    candidates_gff = ARCHAEON.out.candidates_gff

    // Representative model selected from each genomic locus.
    representative_models = ARCHAEON.out.representative_models

    // Extracted representative sequences.
    representative_cds = ARCHAEON.out.representative_cds
    representative_proteins = ARCHAEON.out.representative_proteins

    // Best-reference receptor classifications.
    classified_receptors = ARCHAEON.out.classified_receptors
}


/*
    Define files published from the entry workflow.
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