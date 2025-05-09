#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (params.workflow_opt == 'SYN') {

    ch_input = Channel.fromPath(params.sample_sheet) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.id, file(row.fasta), file(row.gff)) }

    }


include { SYN as SYN } from './workflows/SYN.nf'
include { MULTIQC as MULTIQC } from './workflows/MULTIQC.nf'


workflow {


    if (params.workflow_opt == 'SYN') {

        SYN(ch_input)

        }

    if (params.workflow_opt == 'MULTIQC') {

        MULTIQC()

        }


}