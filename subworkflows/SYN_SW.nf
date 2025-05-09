/*
Subworkflow for assembly QC using Kraken2, minimap2, and blobtools2
Requries set params:


*/

include { BLASTP_ALL as BLASTP_ALL } from '../modules/BLASTP.nf'
include { FASTA_FILT as FASTA_FILT } from '../modules/BIN_SCRIPTS.nf'
include { RENAME_CHR as RENAME_CHR } from '../modules/BIN_SCRIPTS.nf'
include { AGAT_STD as AGAT_STD } from '../modules/AGAT.nf'
include { AGAT_PROT as AGAT_PROT } from '../modules/AGAT.nf'
include { COMBINE_BED as COMBINE_BED } from '../modules/BIN_SCRIPTS.nf'
include { AGAT_GFF2BED as AGAT_GFF2BED } from '../modules/AGAT.nf'


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

        //Combe Protein Fasta & Self-Blast
        BLASTP_ALL(AGAT_PROT.out.prots_only.collect())

        

        //Combine GFFs into BED

        //GFF to BED
        AGAT_GFF2BED(AGAT_STD.out.gff)

        //BED Combine
        COMBINE_BED(AGAT_GFF2BED.out.bed_only.collect())


        //Run McScanX


}