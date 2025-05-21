process AGAT_STD {
    label 'low'
	container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(fa), file(gff)

    output:
	   tuple val(id), path(fa), path("${id}_annot_std.gff3"), emit: gff


    script:

    """
    agat_convert_sp_gxf2gxf.pl --gff ${gff} --output ${id}_annot_std.gff3
    """
}

process AGAT_PROT {
    label 'low'
    container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), path("${id}_proteins.fasta"), emit: prots
       path("${id}_proteins.fasta"), emit: prots_only
       tuple val(id), path(gff), path("${id}_proteins.fasta"), emit: all 


    script:

    """
    agat_sp_extract_sequences.pl \
      --gff ${gff} \
      --fasta ${fa} \
      --type CDS \
      --protein \
      --output ${id}_proteins.fasta

    """
}

process AGAT_GFF2BED {
    label 'low'
    container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), path(fa), path("${id}_annot.bed"), emit: bed
       path("${id}_annot.bed"), emit: bed_only
       tuple val(id), path("${id}_annot.bed"), emit: for_mcscanx


    script:

    """
    agat_convert_sp_gff2bed.pl --gff ${gff} --out ${id}_annot.bed
    """
}

process AGAT_LONGEST_PROT {
    label 'low'
    container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), path("${id}_longest_proteins.fasta"), emit: prots
       path("${id}_longest_proteins.fasta"), emit: prots_only


    script:

    """
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${id}_longest_isoform.gff

    agat_sp_extract_sequences.pl \
      --gff ${id}_longest_isoform.gff \
      --fasta ${fa} \
      --type CDS \
      --protein \
      --output ${id}_longest_proteins.fasta

    """
}

process AGAT_GFF2BED_PAIR {
    label 'low'
    container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(gff1), file(gff2), file(blast)

    output:
       tuple val(id), path("${id}_combined_format.bed"), path(blast), emit: for_mcscanx


    script:

    """
    agat_convert_sp_gff2bed.pl --gff ${gff1} --out gff1_annot.bed
    agat_convert_sp_gff2bed.pl --gff ${gff2} --out gff2_annot.bed

    cat gff1_annot.bed gff2_annot.bed > ${id}_combined.bed
    sort_and_filter_bed.sh ${id}_combined.bed ${id}_combined_format.bed
    """
}