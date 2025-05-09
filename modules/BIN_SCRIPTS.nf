process FASTA_FILT {
    label 'verylow'
	container 'ubuntu:22.04'

    input:
        tuple val(id), file(fa), file(gff)

    output:
	   tuple val(id), file("filt_${fa}""), file(gff), emit: filt_fasta


    script:

    """
    fasta_len_filter.sh ${fa} ${params.fasta_minlen} filt_${fa}
    """
}
process RENAME_CHR {
    label 'verylow'
    container 'ubuntu:22.04'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), file("${id}_genome.fasta"), file("${id}_annot.gff3"), emit: remap


    script:

    """
    generate_contig_map.sh ${fa} contig_map.tsv ${id}
    rename_fasta_headers.sh contig_map.tsv ${fa} ${id}_genome.fasta
    rename_gff3_contigs.sh contig_map.tsv ${gff} ${id}_annot.gff3
    """
}
process COMBINE_BED {
    label 'verylow'
    container 'ubuntu:22.04'

    input:
        file(gffs)

    output:
       file(combined_format.bed), emit: combo_bed


    script:

    """
    cat ${gffs} > combined.bed
    sort_and_filter_bed.sh combined.bed combined_format.bed
    """
}