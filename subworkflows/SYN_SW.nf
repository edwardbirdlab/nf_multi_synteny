/*
Subworkflow for assembly QC using Kraken2, minimap2, and blobtools2
Requries set params:


*/

include { DIAMOND_ALL as DIAMOND_ALL } from '../modules/BLASTP.nf'
include { FASTA_FILT as FASTA_FILT } from '../modules/BIN_SCRIPTS.nf'
include { RENAME_CHR as RENAME_CHR } from '../modules/BIN_SCRIPTS.nf'
include { AGAT_STD as AGAT_STD } from '../modules/AGAT.nf'
include { AGAT_PROT as AGAT_PROT } from '../modules/AGAT.nf'
include { COMBINE_BED as COMBINE_BED } from '../modules/BIN_SCRIPTS.nf'
include { AGAT_GFF2BED as AGAT_GFF2BED } from '../modules/AGAT.nf'
include { MCSCANX as MCSCANX } from '../modules/MCSCANX.nf'
include { AGAT_LONGEST_PROT as AGAT_LONGEST_PROT } from '../modules/AGAT.nf'
include { BUSCO_DB as BUSCO_DB } from '../modules/BUSCO.nf'
include { BUSCO as BUSCO } from '../modules/BUSCO.nf'
include { QUAST as QUAST } from '../modules/QUAST.nf'

workflow SYN_SW {
    take:
        
        input    //    channel: [ val(id), fasta, gff]


    main:

        //Fix fasta files and GFFs filter small contigs


        //Filter small contigs
        FASTA_FILT(input)


        //Rename fasta and gff
        RENAME_CHR(FASTA_FILT.out.filt_fasta)




        //Standardize GFFs

        //AGAT GFF Standard
        AGAT_STD(RENAME_CHR.out.remap)

        //GFFRead Extract Proteins
        AGAT_PROT(AGAT_STD.out.gff)




        //Blast proteins

        //Create pariwise protein set
        protein_ch = AGAT_PROT.out.prots_only
        protein_ch.view{ "Input protein files: $it" }

        pairwise_ch = protein_ch
            .toList()
            .map { files ->
                def labeled = files.indexed().collect { idx, f -> tuple("S${idx+1}", f) }
                labeled.collectMany { t1 ->
                    labeled.collect { t2 ->
                        def label = "${t1[0]}_vs_${t2[0]}"
                        tuple(label, t1[1], t2[1])
                    }
                }
            }

        pairwise_ch.view{ "Pairwise combination: $it" }

        //Running pairwise blasts
        DIAMOND_ALL(pairwise_ch)

        //Collect all blasts
        ch_concatenated_blast = DIAMOND_ALL.out.result
            .collectFile(name: 'all_blast_results.txt', newLine: true)
            .view { file -> "All BLAST results concatenated into: ${file.name}" }

        //Combine GFFs into BED

        //GFF to BED
        AGAT_GFF2BED(AGAT_STD.out.gff)

        //BED Combine
        COMBINE_BED(AGAT_GFF2BED.out.bed_only.collect())


        //Run McScanX
        MCSCANX(ch_concatenated_blast, COMBINE_BED.out.combo_bed)


        //BUSCO

        //Get longest isoform prot seqs
        //AGAT_LONGEST_PROT(AGAT_STD.out.gff)

        //Get Busco DB
        //BUSCO_DB()

        //Run Busco
        //BUSCO(AGAT_LONGEST_PROT.out.prots,BUSCO_DB.out.busco_db)

        //Running Quast
        //QUAST(input)

}