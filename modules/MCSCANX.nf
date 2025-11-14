process MCSCANX {
    label 'mcscanx'
	container 'quay.io/biocontainers/mcscanx:0.1--h9948957_0'

    input:
        tuple val(id), file(gff), file(blast)

    output:
	   tuple path("${id}_mcscanx.html"), path("${id}_combined_input.gff"), path("${id}_combined_input.blast"), path("${id}.collinearity"), emit: output


    script:

    """
    mv ${blast} combined_input_tmp.blast
    awk '{\$1=\$1}1' OFS='\t' combined_input_tmp.blast > combined_input.blast
    mv ${gff} combined_input.gff
    MCScanX combined_input \
        -b ${params.block_pattern} \
        -s ${params.match_size} \
        -e ${params.mcscanx_e_value} \
        -k ${params.match_score} \
        -g ${params.gap_penalty} \
        -m ${params.max_gaps} \
        -w ${params.window_overlap}
    mv *.html ${id}_mcscanx.html
    mv combined_input.gff ${id}_combined_input.gff
    mv combined_input.blast ${id}_combined_input.blast
    mv *.collinearity ${id}.collinearity
    """
}

process MCSCANX_PLEX {
    label 'mcscanx_plex'
    container 'quay.io/biocontainers/mcscanx:0.1--h9948957_0'

    input:
        tuple val(id), file(gff), file(blast)

    output:
       tuple path("${id}_mcscanx.html"), path("${id}_combined_input.gff"), path("${id}_combined_input.blast"), path("${id}.collinearity"), emit: output


    script:

    """
    mv ${blast} combined_input.blast
    mv ${gff} combined_input.gff
    MCScanX combined_input \
        -b ${params.block_pattern} \
        -s ${params.match_size} \
        -e ${params.mcscanx_e_value} \
        -k ${params.match_score} \
        -g ${params.gap_penalty} \
        -m ${params.max_gaps} \
        -w ${params.window_overlap}
    mv *.html ${id}_mcscanx.html
    mv combined_input.gff ${id}_combined_input.gff
    mv combined_input.blast ${id}_combined_input.blast
    mv *.collinearity ${id}.collinearity
    """
}