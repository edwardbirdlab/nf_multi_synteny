include { SYN_SW as SYN_SW } from '../subworkflows/SYN_SW.nf'


workflow SYN {
    take:
        
        genomes    //    channel: [ val(short_name), fasta, gff]

    main:
        SYN_SW(genomes)


}