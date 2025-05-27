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
include { ORTHOFINDER_BG as ORTHOFINDER_BG } from '../modules/ORTHOFINDER.nf'
include { DIAMOND_OF_DB as DIAMOND_OF_DB } from '../modules/BLASTP.nf'
include { DIAMOND_OF as DIAMOND_OF } from '../modules/BLASTP.nf'
include { ORTHOFINDER_BG_RERUN as ORTHOFINDER_BG_RERUN } from '../modules/ORTHOFINDER.nf'
include { COMBINE_BLAST as COMBINE_BLAST } from '../modules/BIN_SCRIPTS.nf'
include { BLAST_RENAME as BLAST_RENAME } from '../modules/BIN_SCRIPTS.nf'
include { AGAT_STATS as AGAT_STATS } from '../modules/AGAT.nf'
include { COMBINE_BED_DUP as COMBINE_BED_DUP } from '../modules/BIN_SCRIPTS.nf'
include { MCSCANX_PLEX as MCSCANX_PLEX } from '../modules/MCSCANX.nf'
include { DIAMOND_PAIR as DIAMOND_PAIR } from '../modules/BLASTP.nf'
include { AGAT_GFF2BED_PAIR as AGAT_GFF2BED_PAIR } from '../modules/AGAT.nf'


workflow SYN_SW {
    take:
        
        input    //    channel: [ val(id), fasta, gff]


    main:

        //Fix fasta files and GFFs filter small contigs


        //Filter small contigs
        // FASTA_FILT(input)


        //Rename fasta and gff
        RENAME_CHR(input)




        //Standardize GFFs

        //AGAT GFF Standard
        AGAT_STD(RENAME_CHR.out.remap)

        //AGAT get stats
        AGAT_STATS(AGAT_STD.out.gff)

        //GFFRead Extract Proteins
        AGAT_LONGEST_PROT(AGAT_STD.out.gff)

        //Combine GFFs into BED

        //GFF to BED
        AGAT_GFF2BED(AGAT_STD.out.gff)

        //BED Combine
        COMBINE_BED(AGAT_GFF2BED.out.bed_only.collect())


        //BUSCO

        //Get longest isoform prot seqs
        AGAT_LONGEST_PROT(AGAT_STD.out.gff)

        //Get Busco DB
        BUSCO_DB()

        //Run Busco
        BUSCO(AGAT_LONGEST_PROT.out.prots,BUSCO_DB.out.busco_db)

        //Running Quast
        QUAST(input)

        //Running Orthofinder3
        ORTHOFINDER_BG(AGAT_LONGEST_PROT.out.prots_only.collect())

        //making diamond dbs
        DIAMOND_OF_DB(ORTHOFINDER_BG.out.fasta.flatten())

        //Channel For Blasts
        ch_diamond_db = DIAMOND_OF_DB.out.db
        ch_diamond_query = ORTHOFINDER_BG.out.fasta.flatten()

        ch_diamond_all = ch_diamond_query.combine(ch_diamond_db)

        //Running Diamond
        DIAMOND_OF(ch_diamond_all)

        //Running orthofinder
        ORTHOFINDER_BG_RERUN(ORTHOFINDER_BG.out.output, DIAMOND_OF.out.result.collect()


        //Making new channel format for McScanX processing
        ch_pairwise = AGAT_LONGEST_PROT.out.all
            .toList()
            .map { list ->
                def combos = []
                for (int i = 0; i < list.size(); i++) {
                    for (int j = i + 1; j < list.size(); j++) {
                        combos << [list[i], list[j]]
                    }
                }
                return combos
            }
            .flatten()
            .buffer(size: 6)
            .filter { i1, g1, p1, i2, g2, p2 -> i1 != i2 }

        DIAMOND_PAIR(ch_pairwise)

        AGAT_GFF2BED_PAIR(DIAMOND_PAIR.out.result)

        MCSCANX_PLEX(AGAT_GFF2BED_PAIR.out.for_mcscanx)

}