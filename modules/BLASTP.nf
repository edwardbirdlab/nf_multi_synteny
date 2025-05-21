process BLASTP_DB {
    label 'midmem'
	container 'ncbi/blast:2.16.0'

    input:
        file(fa)

    output:
	   path("./blast_db_nf"), emit: db


    script:

    """
    mkdir blast_db_nf
    cat ${fa} > all_prots.fasta
    makeblastdb -in all_prots.fasta -dbtype prot -out ./prot_db
    mv prot_db* blast_db_nf
    """
}

process BLASTP {
    label 'low'
    container 'ncbi/blast:2.16.0'

    input:
        tuple file(fa), file(db)

    output:
       path("all_prots_vs_all_prots.blast"), emit: result


    script:

    """
    blastp -query ${fa} \
       -db /${db}/prot_db \
       -out all_prots_vs_all_prots.blast \
       -evalue 1e-10 \
       -max_target_seqs 5 \
       -num_threads ${task.cpus} \
       -outfmt 6
    """
}

process BLASTP_ALL {
    label 'low'
    container 'ncbi/blast:2.16.0'

    input:
        file(fa)

    output:
       path("all_prots_vs_all_prots.blast"), emit: result


    script:

    """
    mkdir blast_db_nf
    cat ${fa} > all_prots.fasta
    makeblastdb -in all_prots.fasta -dbtype prot -out ./prot_db
    mv prot_db* blast_db_nf

    blastp -query all_prots.fasta \
       -db blast_db_nf/prot_db \
       -out all_prots_vs_all_prots.blast \
       -evalue 1e-10 \
       -max_target_seqs 5 \
       -num_threads ${task.cpus} \
       -outfmt 6
    """
}

process DIAMOND_ALL {
    label 'blast'
    container 'buchfink/diamond:version2.1.11'

    input:
        tuple val(id), file(fa1), file(fa2)

    output:
       path("${id}.blast"), emit: result


    script:

    """
    diamond makedb --in ${fa2} -d prots_db

    diamond blastp \
      -q ${fa1} \
      -d prots_db \
      -e ${params.diamond_e_value} \
      --max-target-seqs 5 \
      --outfmt 6 \
      -o ${id}.blast \
      --threads ${task.cpus}
    """
}

process DIAMOND_PAIR {
    label 'blast'
    container 'buchfink/diamond:version2.1.11'

    input:
        tuple val(id1), file(gff1), file(fa1), val(id2), file(gff2), file(fa2)

    output:
       tuple val(sample_combo), path("${sample_combo}_combined_format.bed"), path("${sample_combo}.blast"), emit: result


    script:

    sample_combo = "${id1}_${id2}"

    """
    diamond makedb --in ${fa1} -d prots_db_1
    diamond makedb --in ${fa2} -d prots_db_2

    diamond blastp \
      -q ${fa1} \
      -d prots_db_2 \
      -e ${params.diamond_e_value} \
      --max-target-seqs 5 \
      --outfmt 6 \
      -o ${id1}.blast \
      --threads ${task.cpus}

    diamond blastp \
      -q ${fa2} \
      -d prots_db_1 \
      -e ${params.diamond_e_value} \
      --max-target-seqs 5 \
      --outfmt 6 \
      -o ${id2}.blast \
      --threads ${task.cpus}

    cat ${id1}.blast ${id2}.blast > ${sample_combo}.blast

    cat ${gff1} ${gff2} > ${sample_combo}_combined.bed
    sort_and_filter_bed.sh ${sample_combo}_combined.bed ${sample_combo}_combined_format.bed
    """
}