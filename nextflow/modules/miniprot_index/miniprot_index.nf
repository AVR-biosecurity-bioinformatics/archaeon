process MINIPROT_INDEX {
    tag 'genome'
    conda "${moduleDir}/environment.yml"
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