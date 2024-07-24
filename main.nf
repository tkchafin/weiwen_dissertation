#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Subworkflows 

include { SIMULATE } from './subworkflows/simulate'

workflow {
    // Subworkflow: Simulate sequences
    SIMULATE( params )

    // Subworkflow: 
}