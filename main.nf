#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Subworkflows and module imports
include { SIMULATE } from './subworkflows/simulate'
include { BARCODE_GENERATION } from './subworkflows/barcode_generation'

// Main workflow
workflow {
    // Subworkflow: Simulate sequences
    SIMULATE( params )

    // Subworkflow: Generate barcodes
    BARCODE_GENERATION( params, SIMULATE.out )
}
