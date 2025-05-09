/*
~~~~~~~~~~~~~~~~~~~~~~
Importing subworkflows
~~~~~~~~~~~~~~~~~~~~~~
*/

include { SYN_SN as SYN_SN } from '../subworkflows/SYN_SN.nf'


workflow SYN {
    take:
        
        genomes    //    channel: [ val(short_name), fasta, gff]

    main:
        SYN_SN(genomes)


}