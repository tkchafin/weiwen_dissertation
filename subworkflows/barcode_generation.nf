#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process faToVcf {
    publishDir "${params.outdir}/vcf", mode: 'copy'

    container "pathogengenomics/usher"

    input:
    tuple val(meta), path(align)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    script:
    def prefix = "${meta.id}"
    """
    vcf_file="${prefix}.vcf"
    faToVcf ${align} \${vcf_file}
    """
}

process vcfToBarcodes {
    publishDir "${params.outdir}/barcodes", mode: 'copy'

    container 'tkchafin/pysam:latest'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv"), emit: barcodes

    script:
    def prefix = "${meta.id}"
    """
    vcf_to_barcode.py ${vcf} ${prefix}.tsv --force_unique
    """
}

workflow BARCODE_GENERATION {
    take:
        params
        fasta

    main:
        faToVcf(fasta)
        | vcfToBarcodes
    
    emit:
        vcfToBarcodes.out.barcodes
}
