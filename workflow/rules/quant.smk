# Author: Moe Zardbani
# Quantify reads based on the filtered results from SQANTI3 Filter output
# convert_sq3_output: Modify SQNATI3 Filter classification outptu to make it compatible as 
# input for make_seurat requires additional columns in the classification file.
# make_seurat: Quantify using pigeon PacBio quantitation toolkit.

# Make SQANTI3 filtered classification output compatible with make-seurat
rule convert_sq3_output:
    input: 
        sq3_gtf = "results/sqanti3_filter/{id}/{id}.filtered.gtf",
        pigeon_classify = "results/classify/{id}_classification.txt"
    output: "results/sqanti3_filter/{id}/{id}_RulesFilter_result_classification_seurat.txt"
    conda: "../envs/python_env.yaml"
    benchmark: "results/benchmarks/{id}_convert_sq3_output.benchmark"
    log: "results/logs/{id}_convert_sq3_output.log"
    params:
        sq3_classify = "results/sqanti3_filter/{id}/{id}_RulesFilter_result_classification.txt"
    shell:
        '''
        workflow/scripts/sq3_filter_out_convert.py {params.sq3_classify} \
        {input.pigeon_classify} \
        {output} 2> {log}
        '''

# Quantify SQANTI3 filtered output
rule make_seurat:
    input: "results/sqanti3_filter/{id}/{id}_RulesFilter_result_classification_seurat.txt"
    output: "results/seurat/{id}_sq3/{id}.info.csv"
    conda: "../envs/pbpigeon_env.yaml"
    benchmark: "results/benchmarks/{id}_make_seurat.benchmark"
    log: "results/logs/{id}_make_seurat.log"
    threads: 32
    params:
        dedup = "results/dedup/{id}.fltncc.sorted.dedup.fasta" if config['is_pacbio'] else "results/fasta/{id}.dedup.mapped.miser.corrected.modified.sorted.ont.fasta",
        group = "results/collapse/{id}.collapsed.group.txt",
        out_prefix = "{id}",
        out_dir = "results/seurat/{id}_sq3",
    shell:
        '''
        pigeon make-seurat -j {threads} \
        --log-file {log} \
        --log-level TRACE \
        --keep-ribo-mito-genes \
        --dedup {params.dedup} \
        --group {params.group} \
        -o {params.out_prefix} \
        -d {params.out_dir} \
        {input}
        '''

rule compress_quant:
    input: "results/seurat/{id}_sq3/{id}.info.csv"
    output: "results/seurat/{id}_sq3/genes_seurat/barcodes.tsv.gz"
    conda: "../envs/samtools_env.yaml"
    log: "results/logs/{id}_compress_quant.log"
    benchmark: "results/benchmarks/{id}_compress_quant.benchmark"
    params:
        genes = "results/seurat/{id}_sq3/genes_seurat/genes.tsv",
        features = "results/seurat/{id}_sq3/genes_seurat/features.tsv",
        zip = "results/seurat/{id}_sq3/genes_seurat/*"
    shell:
        '''
        mv {params.genes} {params.features};
        bgzip {params.zip} 2> {log}
        '''