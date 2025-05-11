process ORTHOFINDER_BG {
    label 'lowmem'
	container 'quay.io/biocontainers/orthofinder:3.0.1b1--hdfd78af_0'

    input:
        path(fas, stageAs: 'fastas/')

    output:
        path("OrthoFinder"), emit: output


    script:

    """
    for f in fastas/*_longest_proteins.fasta; do
      newname=\$(basename "\$f" _longest_proteins.fasta).fasta
      cp "\$f" "\$newname"
    done


    mkdir temp_pickle

   orthofinder \\
        -t ${task.cpus} \\
        -a ${task.cpus} \\
        -p temp_pickle \\
        -f fastas \\
        -n ${params.project_name} \\
        -op
    """
}