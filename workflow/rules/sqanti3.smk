## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 4 - SQANTI3

# Workflow description:
# gff2gtf: convert gff to gtf format for SQANTI3 pipeline compatibility
# sqanti3_qc: perform qc and classify transcripts
# sqanti3_filter_rules: filter transcripts based on the rule-based method
# transfrom outputs to prepare for make-seurat with sq3_filter_out_convert.py
# quantify above output with make-seurat

# Prepare input for Sqanti3 QC - Convert collapsed GFF to GTF
rule gff2gtf:
    input: "results/collapse/{id}.collapsed.gff"
    output: "results/collapse/{id}.collapsed.gtf"
    conda: "../envs/gffread_env.yaml"
    benchmark: "results/benchmarks/{id}_gff2gtf.benchmark"
    log: "results/logs/{id}_gff2gtf.log"
    shell:
        '''
        gffread {input} -T -o {output} 2> {log}
        '''

# Generate SQ3 compatible normalised counts from collapsed abundance output
rule prepare_count:
    input: "results/collapse/{id}.collapsed.gtf"
    output: "results/collapse/{id}.collapsed.abundance.sq3.txt"
    conda: "../envs/python_env.yaml"
    benchmark: "results/benchmarks/{id}_prepare_count.benchmark"
    log: "results/logs/{id}_prepare_count.log"
    params:
        input = "results/collapse/{id}.collapsed.abundance.txt",
        wld = "{id}",
    shell:
        '''
        workflow/scripts/gen_sq3_abund.py {params.input} {params.wld} {output} 2> {log}
        '''

# Perform classification and quality control on collapsed transcripts
rule sqanti3_qc:
    input:
        gtf = "results/collapse/{id}.collapsed.gtf",
        fl_count = "results/collapse/{id}.collapsed.abundance.sq3.txt"
    output:
        classification = "results/sqanti3_qc/{id}/{id}_classification.txt"
    conda: "../envs/sq3_env.yaml"
    benchmark: "results/benchmarks/{id}_sqanti3_qc.benchmark"
    log: "results/logs/{id}_sqanti3_qc.log"
    threads: 32
    params:
        cage_peak = config['tss'],
        polyA_motif = config['polyA_motif'],
        # polyA_sites = config['polyA_sites'],
        d = "results/sqanti3_qc/{id}",
        o = "{id}",
        ref_anno = config['anno_gtf'],
        ref_seq = config['ref_fa']
    shell:
        '''
        python workflow/scripts/SQANTI3-5.2.2/sqanti3_qc.py -t {threads} \
        -n 14 \
        --CAGE_peak {params.cage_peak} \
        --polyA_motif_list {params.polyA_motif} \
        --fl_count {input.fl_count} \
        -d {params.d} \
        -o {params.o} \
        --report pdf \
        {input.gtf} \
        {params.ref_anno} \
        {params.ref_seq} 2> {log}
        '''
        # --polyA_peak {params.polyA_sites} \

rule sqanti3_filter_rules:
    input: "results/sqanti3_qc/{id}/{id}_classification.txt"
    output: "results/sqanti3_filter/{id}/{id}.filtered.gtf"
    conda: "../envs/sq3_env.yaml"
    benchmark: "results/benchmarks/{id}_sqanti3_filter_rules.benchmark"
    log: "results/logs/{id}_sqanti3_filter_rules.log"
    params:
        gtf = "results/sqanti3_qc/{id}/{id}_corrected.gtf",
        json = "resources/sq3_filters/filter_exon_level.json",
        d = "results/sqanti3_filter/{id}",
        o = "{id}",
    shell:
        '''
        python workflow/scripts/SQANTI3-5.2.2/sqanti3_filter.py rules \
        --gtf {params.gtf} \
        -j {params.json} \
        -d {params.d} \
        -o {params.o} \
        --skip_report \
        {input} 2> {log}
        '''

# Generate a report on filtering step
rule gen_filter_report:
    input: "results/sqanti3_filter/{id}/{id}.filtered.gtf"
    output: "results/sqanti3_filter/{id}/{id}_SQANTI3_filter_report.pdf"
    conda: "../envs/r_env.yaml"
    benchmark: "results/benchmarks/{id}_gen_filter_report.benchmark"
    log: "results/logs/{id}_gen_filter_report.log"
    params:
        filter_dir = "results/sqanti3_filter/{id}/",
        wld = "{id}",
    shell:
        '''
        Rscript workflow/scripts/SQANTI3_filter_report_moe.R -d {params.filter_dir} \
        -o {params.wld} \
        -u workflow/scripts/SQANTI3-5.2.2/utilities/ \
        -f rules 2> {log}
        '''

# Generate a reads list to use for filtering mapped bam
rule gen_reads_list:
    input: "results/sqanti3_filter/{id}/{id}.filtered.gtf"
    output: "results/pbmm2/{id}_read_names.txt"
    conda: "../envs/python_env.yaml"
    benchmark: "results/benchmarks/{id}_gen_reads_list.benchmark"
    log: "results/logs/{id}_gen_reads_list.log"
    params:
        group = "results/collapse/{id}.collapsed.group.txt",
        inclusion = "results/sqanti3_filter/{id}/{id}_inclusion-list.txt",
        wld = "{id}",
    shell:
        '''
        workflow/scripts/get_read_names.py {params.group} \
        {params.inclusion} \
        {params.wld} \
        {output} 2> {log}
        '''

# Extract filtered reads from mapped bam cell file (filter mapped bam file based on filtered reads)
rule extract_filtered_reads:
    input: "results/pbmm2/{id}_read_names.txt"
    output: "results/pbmm2/{id}.dedup.mapped.miser.corrected.filtered.bam"
    conda: "../envs/samtools_env.yaml"
    benchmark: "results/benchmarks/{id}_extract_filtered_reads.benchmark"
    log: "results/logs/{id}_extract_filtered_reads.log"
    threads: 16
    params:
        bam = "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.bam",
    shell:
        '''
        samtools view -@ {threads} -h -b -N {input} {params.bam} > {output} 2> {log};
        samtools index -@ {threads} {output}
        '''

# Merge filtered transcript models with the reference annotation
rule gffcompare:
    input: "results/sqanti3_filter/{id}/{id}.filtered.gtf"
    output: "results/gffcompare/{id}.annotated.gtf"
    conda: "../envs/gffcompare_env.yaml"
    benchmark: "results/benchmarks/{id}_gffcompare.benchmark"
    log: "results/logs/{id}_gffcompare.log"
    params:
        sorted_anno = "results/sqanti3_filter/{id}/{id}.filtered.sorted.gtf",
        zipped_anno = "results/sqanti3_filter/{id}/{id}.filtered.sorted.gtf.gz",
        ref_anno = config['anno_gtf'],
        out_prefix = "results/gffcompare/{id}"
    shell:
        '''
        sort -k1,1 -k4,4n {input} > {params.sorted_anno};
        gffcompare -r {params.ref_anno} -o {params.out_prefix} {params.sorted_anno} 2> {log};
        bgzip {params.sorted_anno};
        tabix {params.zipped_anno}
        '''
