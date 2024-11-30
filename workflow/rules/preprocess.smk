## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 1

# Workflow description:
# skera: perform segmentation on MAS-Seq reads using MAS-Seq primers list
# lima: remove 10x primers from segmented reads (barcoded samples require barcodes.fa and barcode orientation added to parameters [lima_FAQ](https://lima.how/barcode-design.html))
# tag: clip UMI and barcodes to use for deduplication a later stage (barcode design: https://isoseq.how/umi/umi-barcode-design.html)
# refine: clip polyA tail if persent in sample

# Perform segmentation on MAS-Seq reads if reads are NOT already segmented
if not config['is_segmented']:
    rule skera:
        input: lambda wildcards: input_files[wildcards.id]
        output: "results/skera/{id}.segmented.bam"
        conda: "../envs/pbskera_env.yaml"
        benchmark: "results/benchmarks/{id}_skera.benchmark"
        threads: 16
        params:
            log = "results/logs/{id}_skera.log",
            log_level = "TRACE",
            mas_primers = config['mas_primers']
        shell:
            '''
            skera split -j {threads} \
            --log-level {params.log_level} \
            --log-file {params.log} \
            {input} {params.mas_primers} {output}
            '''

# Remove 10x primers and re-orient reads to 5'->3'
rule lima:
    input: lambda wildcards: input_files[wildcards.id] if config['is_segmented'] else "results/skera/{id}.segmented.bam"
    output: "results/lima/{id}.fl.5p--3p.bam"
    conda: "../envs/lima_env.yaml"
    benchmark: "results/benchmarks/{id}_lima.benchmark"
    threads: 16
    params: 
        log = "results/logs/{id}_lima.log",
        log_level = "TRACE",
        primers = config['primers'],
        outfile = "results/lima/{id}.fl.bam"
    shell:
        '''
        lima --isoseq {input} \
        {params.primers} \
        {params.outfile} \
        --per-read \
        --peek-guess \
        --log-level {params.log_level} \
        --log-file {params.log} \
        -j {threads}
        '''

# Clip UMI/barcodes and store for downstream deduplication (barcode design should be specified: https://isoseq.how/umi/umi-barcode-design.html e.g. 10X 3' Chromium kit design = T-12U-16B)
rule tag:
    input: "results/lima/{id}.fl.5p--3p.bam"
    output: "results/tag/{id}.flt.bam"
    conda: "../envs/isoseq_env.yaml"
    benchmark: "results/benchmarks/{id}_tag.benchmark"
    threads: 16
    params:
        log = "results/logs/{id}_tag.log",
        log_level = "TRACE",
        # Adjust design based on experiment
        design = config['lib_design']
    shell:
        '''
        isoseq tag --design {params.design} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {threads} \
        {input} \
        {output}
        '''
# rm {input.lima}

# Trim polyA and remove concatenated reads
# ***Note that output of this step is an optional input for `collapse` (bulk ISO-Seq data)
rule refine:
    input: "results/tag/{id}.flt.bam"
    output: "results/refine/{id}.fltnc.bam"
    conda: "../envs/isoseq_env.yaml"
    benchmark: "results/benchmarks/{id}_refine.benchmark"
    threads: 16
    params:
        primer = config['primers'],
        log = "results/logs/{id}_refine.log",
        log_level = "TRACE"
    shell:
        '''
        isoseq refine \
        --require-polya \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {threads} \
        {input} \
        {params.primer} \
        {output}
        '''
# rm {input.tag}