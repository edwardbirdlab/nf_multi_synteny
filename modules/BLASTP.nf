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
       -db /blast_db_nf/prot_db \
       -out all_prots_vs_all_prots.blast \
       -evalue 1e-10 \
       -max_target_seqs 5 \
       -num_threads ${task.cpus} \
       -outfmt 6
    """
}