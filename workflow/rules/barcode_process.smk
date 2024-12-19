## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: MoMoTPark
# PART 2

# Workflow description:
# correct: perform barcode correction to recover more reads and make deduplication more effective
# bcstats: get cell calling statistics to adjust correction methods if required
# plot_knees: generate cell calling QC knee plots

rule correct:
    input: "results/refine/{id}.fltnc.bam"
    output: "results/barcode_corrected/{id}.fltncc.bam"
    conda: "../envs/isoseq_env.yaml"
    log: "results/logs/{id}_correct.log"
    benchmark: "results/benchmarks/{id}_correct.benchmark"
    threads: 16
    params:
        barcodes = config['barcode_whitelist']
    shell:
        '''
        isoseq correct -j {threads} \
        --method knee \
        --log-file {log} \
        --log-level TRACE \
        --barcodes {params.barcodes} \
        {input} \
        {output}
        '''

# Calculate cell calling statistics
rule bcstats:
    input: "results/barcode_corrected/{id}.fltncc.bam"
    output: "results/bcstats/{id}.bcstats.tsv"
    conda: "../envs/isoseq_env.yaml"
    log: "results/logs/{id}_bcstats.log"
    benchmark: "results/benchmarks/{id}_bcstats.benchmark"
    threads: 16
    params:
        log_level = "TRACE",
        json = "results/bcstats/{id}.bcstats.json"
    shell:
        '''
        isoseq bcstats -j {threads} \
        --method knee \
        --log-file {log} \
        --log-level TRACE \
        -o {output} \
        --json {params.json} \
        {input}
        ''' 

# Generate cell calling QC knee plots
rule plot_knees:
    input: "results/bcstats/{id}.bcstats.tsv"
    output: "results/bcstats/{id}.knee.png"
    conda: "../envs/plotknees_env.yaml"
    log: "results/logs/{id}_plot_knees.log"
    benchmark: "results/benchmarks/{id}_plot_knees.benchmark"
    params:
        percentile = 95,
        output = "results/bcstats/{id}",
    shell:
        '''
        python workflow/scripts/plot_knees.py \
        --tsv {input} \
        --output {params.output} \
        --estimate_percentile {params.percentile} 2> {log}
        '''