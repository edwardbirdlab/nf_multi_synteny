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
       tuple val(sample_combo), path(gff1), path(gff2), path("${sample_combo}.blast"), emit: result
       tuple val(id1), path("${id1}_vs_${id2}.blast"), emit: result_sp1
       tuple val(id2), path("${id2}_vs_${id1}.blast"), emit: result_sp2


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
      -o ${id1}_vs_${id2}.blast \
      --threads ${task.cpus}

    diamond blastp \
      -q ${fa2} \
      -d prots_db_1 \
      -e ${params.diamond_e_value} \
      --max-target-seqs 5 \
      --outfmt 6 \
      -o ${id2}_vs_${id1}.blast \
      --threads ${task.cpus}

    cat ${id1}_vs_${id2}.blast ${id2}_vs_${id1}.blast > ${sample_combo}.blast

    cat ${gff1} ${gff2} > ${sample_combo}_combined.bed
    sort_and_filter_bed.sh ${sample_combo}_combined.bed ${sample_combo}_combined_format.bed
    "

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
        tuple file(fa), file(db)

    output:
       path("*.txt.gz"), emit: result


    script:

    """
    Q_ID=\$(basename ${fa} .fa | grep -o '[0-9]\\+')
    DB_ID=\$(basename ${db} .dmnd | grep -o '[0-9]\\+')
    DB_NAME=\$(basename ${db} .dmnd)

    OUTFILE="Blast\${Q_ID}_\${DB_ID}.txt"
    
    diamond blastp -d \$DB_NAME -q ${fa} -o \$OUTFILE --more-sensitive -p ${task.cpus} --quiet -e 0.001 --compress 1
    """
}