process BUSCO_DB {
    label 'lowmem'
    container 'ezlabgva/busco:v5.8.2_cv1'
        
    output:
        path("./busco_downloads"), emit: busco_db
        path("versions.yml"), emit: versions

    script:

    """
    busco --download ${params.busco_lineage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco -v 2>&1 | sed -e "s/BUSCO //g")
    END_VERSIONS 
    """
}

process BUSCO {
    label 'lowmem'
    container 'ezlabgva/busco:v5.8.2_cv1'

    input:
        tuple val(sample), file(fasta)
        file(busco_db)

    output:
        path("./${sample}_busco"), emit: busco_results
        path("versions.yml"), emit: versions

    script:

    """
    busco -i ${fasta} -o ${sample}_busco -m genome --offline --download_path ${busco_db} --lineage_dataset ${params.busco_lineage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco -v 2>&1 | grep "BUSCO" | sed -e "s/BUSCO //g")
    END_VERSIONS  
    """
}