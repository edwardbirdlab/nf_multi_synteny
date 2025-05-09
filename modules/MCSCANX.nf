process MCSCANX {
    label 'midmem'
	container 'quay.io/biocontainers/mcscanx:0.1--h9948957_0'

    input:
        file(blast)
        file(gff)

    output:
	   path("*.html"), path("*.gff"), path("*.fasta"), path("*.collinearity"), emit: output


    script:

    """
    mv ${blast} combined_input.blast
    mv ${gff} combined_input.gff
    MCScanX combined_input \
        -b ${params.block_pattern} \
        -s ${params.match_size} \
        -e ${params.e_value} \
        -k ${params.match_score} \
        -g ${params.gap_penalty} \
        -m ${params.max_gaps} \
        -w ${params.window_overlap}
    """
}