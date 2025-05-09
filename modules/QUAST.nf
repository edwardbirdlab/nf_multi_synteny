process QUAST {
    label 'lowmem'
    container 'ebird013/quast:5.2.0'

    input:
        tuple val(id), file(fasta), file(gff)

    output:
        path("./${id}"), emit: quast_results
        path("versions.yml"), emit: versions

    script:

    """
    mv ${fasta} ${id}.fasta
    quast.py -o ${id} ${id}.fasta --threads ${task.cpus} --min-contig 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Quast: \$(quast --version 2>&1)
    END_VERSIONS 
    """
}