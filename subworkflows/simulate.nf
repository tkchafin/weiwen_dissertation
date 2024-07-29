#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def generateSdValues(min, max, step) {
    def sd_values = []
    for (float sd = min; sd <= max; sd += step) {
        sd_values << sd
    }
    return sd_values
}

process simulateTrees {
    publishDir "${params.outdir}/trees", mode: 'copy'

    container 'tkchafin/dendropy:latest'

    input:
    val sd_fixed

    output:
    path("*.tre"), emit: tree_files

    script:
    """
    simulate_trees.py --num_tips ${params.num_tips} --sd_fixed ${sd_fixed} --reps ${params.reps} --out ./sim --seed ${params.seed}
    """
}

process selectAndCleanTrees {
    publishDir "${params.outdir}/selected_trees", mode: 'copy'

    container 'tkchafin/dendropy:latest'

    input:
    path tree_files

    output:
    path("selected*.tre"), emit: selected_tree_files
    path("*.pdf"), emit: report_files
    path("*selected_trees.txt"), emit: selected_trees_list

    script:
    """
    select_trees.py --input ./ --out ./selected
    """
}

process rescaleTrees {
    publishDir "${params.outdir}/scaled_trees", mode: 'copy'

    container 'tkchafin/dendropy:latest'

    input:
    tuple path(tree_file), val(scale_factor)

    output:
    path "*scale*.tre", emit: rescaled_tree_files

    script:
    """
    rescale_branch_length.py --input_tree ${tree_file} --scale ${scale_factor} --output_dir ./
    """
}

process generateReferenceGenome {
    publishDir "${params.outdir}/references", mode: 'copy'

    input:
    tuple val(meta), path(rescaled_tree_file)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    script:
    def prefix = "${meta.id}"
    """
    BASENAME=\$(basename ${rescaled_tree_file} .tre)
    REF_FILE=${prefix}.fasta
    echo ">reference" > \${REF_FILE}
    cat /dev/urandom | tr -dc 'ACGT' | fold -w ${params.genome_length} | head -n 1 >> \${REF_FILE}
    echo "${rescaled_tree_file} \${REF_FILE}"
    """
}

process phastSim {
    publishDir "${params.outdir}/sequences", mode: 'copy'

    container 'tkchafin/phastsim:latest'

    input:
    tuple val(meta), path(tree_file), path(reference)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.info"), emit: info

    script:
    def prefix = "${meta.id}"
    """
    phastSim --outpath ./ --reference ${reference} --treeFile ${tree_file} \
        --categoryRates 1.0 1.1 1.5 2.0 --categoryProbs 0.9 0.01 0.04 0.05 \
        --outputFile ${prefix} --createFasta --createInfo
    """
}

process checkAlignments {
    publishDir "${params.outdir}/reports", mode: 'copy'

    container 'tkchafin/numpy:latest'

    input:
    tuple val(meta), path(fasta_file)

    output:
    path("*_report.txt"), emit: report

    script:
    def prefix = "${meta.id}"
    """
    check_alignments.py --fasta ${fasta_file} --tips ${params.num_tips} --prefix ${prefix}

    num_unique=\$(grep "Number of unique sequences:" ${prefix}_report.txt | awk '{print \$5}')
    if [ "\$num_unique" -lt "${params.num_tips}" ]; then
        echo "Warning: Number of unique sequences (\$num_unique) is below the expected number of tips (${params.num_tips})." >> ${prefix}_report.txt
    fi
    """
}

process collectReports {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path(report_files)

    output:
    path("sequence_report.txt")

    script:
    """
    echo "Meta ID, Number of Unique Sequences, Number of Duplicate Sequences, Minimum Hamming Distance, Maximum Hamming Distance, Average Hamming Distance" > sequence_report.txt

    for report in ${report_files}; do
        meta_id=\$(basename \$report _report.txt)
        num_unique=\$(grep "Number of unique sequences:" \$report | awk '{print \$5}')
        num_duplicates=\$(grep "Number of duplicate sequences:" \$report | awk '{print \$5}')
        min_distance=\$(grep "Minimum Hamming distance:" \$report | awk '{print \$4}')
        max_distance=\$(grep "Maximum Hamming distance:" \$report | awk '{print \$4}')
        avg_distance=\$(grep "Average Hamming distance:" \$report | awk '{print \$4}')
        
        echo "\${meta_id}, \${num_unique}, \${num_duplicates}, \${min_distance}, \${max_distance}, \${avg_distance}" >> sequence_report.txt
    done
    """
}

// process makeUnique {
//     publishDir "${params.outdir}/sequences_clean", mode: 'copy'

//     //container 'tkchafin/numpy:latest'

//     input:
//     tuple val(meta), path(fasta_file)

//     output:
//     tuple val(meta), path("*.validated.fasta"), emit: fasta

//     script:
//     def prefix = "${meta.id}"
//     """
//     make_unique_fasta.py --fasta ${fasta_file} --prefix ${prefix}
//     """
// }

workflow SIMULATE {
    take:
        params

    main:
        Channel
            .from(generateSdValues(params.sd_min, params.sd_max, params.sd_step))
            .set { sd_values_channel }

        Channel
            .from(params.scale_factors)
            .set { scale_factors_channel }

        // Simulate trees
        simulateTrees(sd_values_channel)

        // Select trees
        trees = simulateTrees.out.tree_files.flatten().collect()
        selectAndCleanTrees(trees)

        // Combine trees and scale factors
        selectAndCleanTrees.out.selected_tree_files
            | flatMap { it }
            | map { file(it) }
            | combine (scale_factors_channel)
            | map { tree, scale -> tuple(tree, scale) }
            | rescaleTrees


        // Add meta to rescaleTrees output channel
        rescaleTrees.out.rescaled_tree_files
            | map { file -> tuple(id: file.baseName, file) }
            | set { ch_rescaled_trees }

        // Generate random reference genomes
        generateReferenceGenome( ch_rescaled_trees )
        
        // Generate sequences
        phastSim(
            ch_rescaled_trees
                | join( generateReferenceGenome.out.fasta )
        )
        
        // // Check alignments
        // checkAlignments(phastSim.out.fasta)
        //     .flatMap { it }
        //     .collect { it }
        //     .set { report_files }

        // // Collect reports
        // collectReports( report_files )


    emit:
        phastSim.out.fasta
}
