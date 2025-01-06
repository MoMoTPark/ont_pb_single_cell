# Author: Moe
# ONT specific rules
# Nextflow pipeline: https://github.com/epi2me-labs/wf-single-cell

# Run Nextflow pipeline for ONT data
rule ont_pipeline:
    input: lambda wildcards: input_files[wildcards.id]
    output: "results/{id}/tagged.bam"
    log: "results/logs/{id}_ont_pipeline.log"
    benchmark: "results/benchmarks/{id}_ont_pipeline.benchmark"
    threads: 128
    conda: "../envs/nextflow_env.yaml"
    params:
        sample = "{id}",
        outdir = "results/",
        workdir = ".nf_tmp/{id}/",
        expected_cells = lambda wildcards: expected_cells[wildcards.id],
        mito_prefix = config['mito_prefix'],
        ont_ref = config['ont_ref'],
        sc_kit = config['sc_kit']
    shell:
        '''
        nextflow run resources/sc_ont_wf/main.nf \
        -resume \
        -profile singularity \
        -c resources/sc_ont_wf/numba.config \
        -w {params.workdir} \
        --threads {threads} \
        --matrix_min_genes 0 \
        --matrix_min_cells 0 \
        --matrix_max_mito 0 \
        --expected_cells {params.expected_cells} \
        --sample {params.sample} \
        --fastq '{input}' \
        --kit '{params.sc_kit}' \
        --ref_genome_dir '{params.ont_ref}' \
        --out_dir {params.outdir} \
        --full_length_only false \
        --mito_prefix '{params.mito_prefix}' 2> {log}
        '''

# Add RG header and isoseq collapse compatible tags
rule modify_ont_bam:
    # Add only_if: in case expand conditional does not work as expected
    input: "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.bam"
    output: "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.ont.bam"
    conda: "../envs/pysam_env.yaml"
    log: "results/logs/{id}_modify_ont_bam.log"
    benchmark: "results/benchmarks/{id}_modify_ont_bam.benchmark"
    shell:
        '''
        python workflow/scripts/add_new_tag.py {input} {output} 2> {log}
        '''

# Generate FASTA from mapped bam
rule gen_ont_fasta:
    input: "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.ont.bam"
    output: "results/fasta/{id}.dedup.mapped.miser.corrected.modified.sorted.ont.fasta"
    conda: "../envs/samtools_env.yaml"
    log: "results/logs/{id}_gen_ont_fasta.log"
    benchmark: "results/benchmarks/{id}_gen_ont_fasta.benchmark"
    threads: 16
    shell:
        '''
        samtools fasta -@ {threads} \
        {input} > {output} 2> {log}
        '''
