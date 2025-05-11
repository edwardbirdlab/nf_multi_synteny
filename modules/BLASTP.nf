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

process DIAMOND_OF_DB {
    label 'blast'
    container 'buchfink/diamond:version2.1.11'

    input:
        file(fa)

    output:
       path("*.dmnd"), emit: db


    script:

    """
    BASENAME=\$(basename ${fa} .fa)
    
    diamond makedb --in ${fa} -d diamondDB\${BASENAME}
    """
}

process DIAMOND_OF {
    label 'blast'
    container 'buchfink/diamond:version2.1.11'

    input:
        tupel file(fa), file(db)

    output:
       path("*.txt.gz"), emit: result


    script:

    """
    Q_ID=\$(basename ${query} .fa | grep -o '[0-9]\\+')
    DB_ID=\$(basename ${db} .dmnd | grep -o '[0-9]\\+')
    DB_NAME=\$(basename ${db} .dmnd

    OUTFILE="Blast\${Q_ID}_\${DB_ID}.txt"
    
    diamond blastp -d \$DB_NAME -q ${fa} -o \$OUTFILE --more-sensitive -p 1 --quiet -e 0.001 --compress 1
    """
}