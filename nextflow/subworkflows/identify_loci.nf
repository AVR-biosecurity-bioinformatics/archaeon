
/*
    Genotype samples using bcftools
*/


//// import modules
include { MINIPROT_INDEX                  } from '../modules/miniprot_index/miniprot_index' 
include { MINIPROT_ALIGN                  } from '../modules/miniprot_align/miniprot_align' 
include { FILTER_CANDIDATES               } from '../modules/filter_candidates/filter_candidates' 
include { SELECT_REPRESENTATIVE_LOCI      } from '../modules/select_representative_loci/select_representative_loci' 
include { CLASSIFY_REPRESENTATIVE_LOCI    } from '../modules/classify_representative_loci/classify_representative_loci' 
include { GFFREAD_EXTRACT                 } from '../modules/gffread_extract/gffread_extract' 


workflow IDENTIFY_LOCI {

    take:
    genome
    references
    metadata

    main:

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

    emit:

    // Raw Miniprot outputs.
    miniprot_gff = MINIPROT_ALIGN.out.gff
    miniprot_paf = MINIPROT_ALIGN.out.paf
    miniprot_models = MINIPROT_ALIGN.out.models

    // Models passing coverage and identity thresholds.
    candidates_tsv = FILTER_CANDIDATES.out.summary
    candidates_gff = FILTER_CANDIDATES.out.gff

    // One provisional representative model per genomic locus.
    representative_models = SELECT_REPRESENTATIVE_LOCI.out.models

    // CDS and translated protein sequences for representative models.
    representative_cds = GFFREAD_EXTRACT.out.cds
    representative_proteins = GFFREAD_EXTRACT.out.proteins

    // Best-reference classification for each representative locus.
    classified_receptors = CLASSIFY_REPRESENTATIVE_LOCI.out.classifications
}