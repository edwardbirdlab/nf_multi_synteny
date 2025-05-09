/*
~~~~~~~~~~~~~~~~~~~~~~
Importing subworkflows
~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC_SW as MULTIQC_SW } from '../subworkflows/MULTIQC_SW.nf'


workflow MULTIQC {

    main:
        MULTIQC_SW()

}