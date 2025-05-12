process FASTA_FILT {
    label 'verylow'
	container 'ubuntu:22.04'

    input:
        tuple val(id), file(fa), file(gff)

    output:
	   tuple val(id), path("filt_${fa}"), path(gff), emit: filt_fasta


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
       tuple val(id), path("${id}_genome.fasta"), path("${id}_annot.gff3"), emit: remap
       path("${id}_contig_map.tsv"), emit: map


    script:

    """
    generate_contig_map.sh ${fa} ${id}_contig_map.tsv ${id}_
    rename_fasta_headers.sh ${id}_contig_map.tsv ${fa} ${id}_genome.fasta
    rename_gff3_contigs.sh ${id}_contig_map.tsv ${gff} ${id}_annot.gff3
    """
}
process COMBINE_BED {
    label 'verylow'
    container 'ubuntu:22.04'

    input:
        file(gffs)

    output:
       path("combined_format.bed"), emit: combo_bed


    script:

    """
    cat ${gffs} > combined.bed
    sort_and_filter_bed.sh combined.bed combined_format.bed
    """
}
process COMBINE_BLAST {
    label 'verylow'
    container 'ubuntu:22.04'

    input:
        file(blasts)

    output:
       path("combined_input.blast"), emit: combo


    script:

    """
    gunzip *.gz
    cat *.txt > combined_input.blast
    """
}
process PCG_COUNT {
    label 'verylow'
    container 'ubuntu:22.04'

    input:
        tuple val(id), file(prots)

    output:
       path("${id}_prot_count.txt"), emit: counts


    script:

    """
    grep -c '^>' ${prots} > ${id}_prot_count.txt
    """
}