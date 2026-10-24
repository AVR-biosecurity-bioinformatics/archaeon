params {

    // Mandatory inputs
    genome: Path
    references: Path
    slurm_account: String                      // SLURM account for submitted jobs.

    // Optional inputs
    outdir: String = "results"
    min_coverage: Float = 0.5f
    min_identity: Float = 0.3f
    max_intron: Integer = 200000
    max_hits: Integer = 50
    secondary_ratio: Float = 0.5f
    translation_table: Integer = 1

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
    conda 'bioconda::miniprot=0.18'
    cpus 4
    memory '8 GB'

    input:
    path index, stageAs: 'genome.mpi'
    path references, stageAs: 'references.faa'

    output:
    path 'miniprot.gff3', emit: gff
    script:
    """
    miniprot -t ${task.cpus} \
     -j 2 \
     -G ${params.max_intron} \
     -N ${params.max_hits} \
     --outn ${params.max_hits} \
     -p ${params.secondary_ratio} \
     --outs ${params.secondary_ratio} \
     --outc ${params.min_coverage} \
     --gff genome.mpi references.faa \
     > miniprot.gff3
    """
}

process EXTRACT_CANDIDATES {
    tag 'extract and QC'
    conda 'bioconda::python=3.11 bioconda::biopython=1.85'
    cpus 1
    memory '4 GB'

    input:
    path genome, stageAs: 'genome.fa'
    path references, stageAs: 'references.faa'
    path gff, stageAs: 'miniprot.gff3'

    output:
    path 'candidates.tsv', emit: summary
    path 'candidates.gff3', emit: gff
    path 'candidates.cds.fna', emit: cds
    path 'candidates.proteins.faa', emit: proteins
    path 'qc.json', emit: qc
    script:
    """
    python extract_candidates.py \
        --genome genome.fa \
        --references references.faa \
        --gff miniprot.gff3 \
        --min-coverage ${params.min_coverage} \
        --min-identity ${params.min_identity} \
        --translation-table ${params.translation_table}
    """
}

/*
    Nextflow workflow 
    This is where processes are connected and pass their outputs to downstream steps
*/

workflow {

    main:
    genome = channel.value(params.genome)
    references = channel.value(params.references)

    MINIPROT_INDEX(genome)

    MINIPROT_ALIGN(
        MINIPROT_INDEX.out.index,
        references
    )

    EXTRACT_CANDIDATES(
        genome,
        references,
        MINIPROT_ALIGN.out.gff
    )

    publish:

    miniprot_gff = MINIPROT_ALIGN.out.gff

    candidates_tsv = EXTRACT_CANDIDATES.out.summary
    candidates_gff = EXTRACT_CANDIDATES.out.gff
    candidates_cds = EXTRACT_CANDIDATES.out.cds
    candidates_proteins = EXTRACT_CANDIDATES.out.proteins
    candidates_qc = EXTRACT_CANDIDATES.out.qc
}

/*
    Publish outputs 
    This is where published outputs of the workflow are sent to your output directory
*/

output {
    miniprot_gff { path 'alignments' }
    candidates_tsv { path 'candidates' }
    candidates_gff { path 'candidates' }
    candidates_cds { path 'candidates' }
    candidates_proteins { path 'candidates' }
    candidates_qc { path 'candidates' } 
}