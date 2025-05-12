process AGAT_STD {
    label 'low'
	container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(fa), file(gff)

    output:
	   tuple val(id), path(fa), path("${id}_annot_std.gff3"), emit: gff
       file("versions.yml"), emit: versions


    script:

    """
    agat_convert_sp_gxf2gxf.pl --gff ${gff} --output ${id}_annot_std.gff3

      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_statistics.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
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
       file("versions.yml"), emit: versions


    script:

    """
    agat_sp_extract_sequences.pl \
      --gff ${gff} \
      --fasta ${fa} \
      --type CDS \
      --protein \
      --output ${id}_proteins.fasta

      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_statistics.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
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
       file("versions.yml"), emit: versions


    script:

    """
    agat_convert_sp_gff2bed.pl --gff ${gff} --out ${id}_annot.bed

      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_statistics.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
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
       file("versions.yml"), emit: versions


    script:

    """
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${id}_longest_isoform.gff

    agat_sp_extract_sequences.pl \
      --gff ${id}_longest_isoform.gff \
      --fasta ${fa} \
      --type CDS \
      --protein \
      --output ${id}_longest_proteins.fasta

      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_statistics.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}

process AGAT_STATS {
    label 'low'
    container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'

    input:
        tuple val(id), file(fa), file(gff)

    output:
       tuple val(id), path("${id}_agat_stats.fasta"), emit: stats
       file("versions.yml"), emit: versions


    script:

    """
    agat_sp_statistics.pl \
      --gff ${gff} \
      --f ${fa} \
      --output ${id}_agat_stats.fasta


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_statistics.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}