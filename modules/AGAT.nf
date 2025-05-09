process AGAT_STD {
    label 'low'
	container 'ebird013/agat:latest'

    input:
        tuple val(id), file(fa), file(gff)

    output:
	   tuple val(id), file(fa), file(${id}_annot_std.gff3), emit: gff


    script:

    """
    agat_convert_sp_gff2gff3.pl --gff ${gff} > ${id}_annot_std.gff3
    """
}

process AGAT_PROT {
    label 'low'
    container 'ebird013/agat:latest'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), file(${id}_proteins.fasta), emit: prots
       file(${id}_proteins.fasta), emit: prots_only


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
    container 'ebird013/agat:latest'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), file(fa), file(${id}_annot.bed), emit: bed
       file(${id}_annot.bed), emit: bed_only


    script:

    """
    agat_convert_gff2bed.pl --gff ${gff} --bed ${id}_annot.bed
    """
}