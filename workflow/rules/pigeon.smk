## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 4 - Pigeon

# Workflow description:
# pigeon_prepare: prepare input files for pigeon to classify and filter transcripts
# pigeon_classify: classification of transcripts into SQANTI3 defined categories
# pigeon_filter: filter transcripts to reduce false positives
# pigeon_report: prepare a report for classification and filtering steps
# make_seurat: create a barcode-feature matrix compatible with Seurat as input

# Sort collapse output gff
rule pigeon_prepare:
    input: "results/collapse/{id}.collapsed.gff"
    output: "results/collapse/{id}.collapsed.sorted.gff"
    conda: "../envs/pbpigeon_env.yaml"
    benchmark: "results/benchmarks/{id}_pigeon_prepare.benchmark"
    log: "results/logs/{id}_pigeon_prepare.log"
    shell:
        '''
        pigeon prepare --log-file {log} \
        --log-level TRACE \
        {input}
        '''

# Classify transcripts into categories: https://isoseq.how/classification/categories
rule pigeon_classify:
    input: "results/collapse/{id}.collapsed.sorted.gff"
    output: "results/classify/{id}_classification.txt"
    conda: "../envs/pbpigeon_env.yaml"
    benchmark: "results/benchmarks/{id}_pigeon_classify.benchmark"
    log: "results/logs/{id}_pigeon_classify.log"
    threads: 32
    params:
        out_dir = "results/classify",
        prefix = "{id}",
        anno = config['anno_gtf'],
        ref = config['ref_fa'],
        polyA = config['polyA_motif'],
        fl_counts = "results/collapse/{id}.collapsed.abundance.txt",
    shell:
        '''
        pigeon index {params.anno};
        pigeon classify -j {threads} \
        --log-file {log} \
        --log-level TRACE \
        -d {params.out_dir} \
        -o {params.prefix} \
        --poly-a {params.polyA} \
        --flnc {params.fl_counts} \
        {input} \
        {params.anno} \
        {params.ref}
        '''
# --gene-id \