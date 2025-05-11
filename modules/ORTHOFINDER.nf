process ORTHOFINDER_BG {
    label 'lowmem'
	container 'quay.io/biocontainers/orthofinder:3.0.1b1--hdfd78af_0'

    input:
        path(fas, stageAs: 'fastas/')

    output:
        path("OrthoFinder"), emit: output
        path("OrthoFinder/Results_${params.project_name}/WorkingDirectory/*.fa"), emit: fasta
        path("OrthoFinder/Results_${params.project_name}/WorkingDirectory/SequenceIDs.txt"), emit: seqids
        path("OrthoFinder/Results_${params.project_name}/WorkingDirectory/SpeciesIDs.txt"), emit: spids


    script:

    """
    for f in fastas/*_longest_proteins.fasta; do
      newname="fastas/\$(basename "\$f" _longest_proteins.fasta).fasta"
      mv "\$f" "\$newname"
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

process ORTHOFINDER_BG_RERUN {
    label 'orthofinder'
    container 'quay.io/biocontainers/orthofinder:3.0.1b1--hdfd78af_0'

    input:
        path(of_dir)
        path(blast, stageAs: 'OrthoFinder/Results_${params.project_name}/WorkingDirectory')

    output:
        path("OrthoFinder"), emit: output


    script:

    """
    mkdir temp_pickle

   orthofinder \\
        -t ${task.cpus} \\
        -a ${task.cpus} \\
        -p temp_pickle \\
        -b OrthoFinder \\
        -n ${params.project_name} \\
        -op
    """
}