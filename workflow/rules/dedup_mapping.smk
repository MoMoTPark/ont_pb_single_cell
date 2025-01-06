## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 3

# Workflow description:
# sort: sort bam files with corrected barcodes by `CB` tag to prepare input as groupdedup
# dedup: deduplicate reads based on UMIs grouped by corrected cell barcodes
# pbmm2: map reads to indexed reference genome
# miser: correct any micro-exon misalignments

# Sort corrected bam file by CB tag 
rule sort:
    input: "results/barcode_corrected/{id}.fltncc.bam"
    output: "results/sorted/{id}.fltncc.sorted.bam"
    conda: "../envs/samtools_env.yaml"
    log: "results/logs/{id}_sort.log"
    benchmark: "results/benchmarks/{id}_sort.benchmark"
    threads: 16
    shell:
        '''
        samtools sort -@ {threads} \
        -t CB \
        {input} \
        -O BAM \
        -o {output} 2> {log}
        '''

# Deduplicate reads based on UMI grouped by cell barcodes
rule dedup:
    input: "results/sorted/{id}.fltncc.sorted.bam"
    output: "results/dedup/{id}.fltncc.sorted.dedup.bam"
    conda: "../envs/isoseq_env.yaml"
    benchmark: "results/benchmarks/{id}_dedup.benchmark"
    threads: 32
    params:
        log = "results/logs/{id}_dedup.log",
        log_level = "TRACE"
    shell:
        '''
        isoseq groupdedup -j {threads} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        --keep-non-real-cells \
        {input} \
        {output}
        '''

# Map reads to reference with pbmm2
rule pbmm2:
    input: "results/dedup/{id}.fltncc.sorted.dedup.bam"
    output: "results/pbmm2/{id}.dedup.mapped.bam"
    conda: "../envs/pbmm2_env.yaml"
    benchmark: "results/benchmarks/{id}_pbmm2.benchmark"
    threads: 32
    params:
        log = "results/logs/{id}_pbmm2.log",
        log_level = "TRACE",
        ref = config['ref_mmi']
    shell:
        '''
        pbmm2 align --preset ISOSEQ \
        --sort \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {threads} \
        {params.ref} \
        {input} \
        {output}
        '''

# Correct microexon misallignment with MisER
rule miser:
    input: "results/pbmm2/{id}.dedup.mapped.bam" if config['is_pacbio'] else "results/{id}/tagged.bam"
    output: "results/pbmm2/{id}.dedup.mapped.miser.corrected.bam"
    conda: "../envs/miser_env.yaml"
    benchmark: "results/benchmarks/{id}_miser.benchmark"
    log: "results/logs/{id}_miser.log"
    threads: 32
    params:
        out_region = "results/pbmm2/{id}.dedup.mapped.miser.regions.out.txt",
        ref = config['ref_fa'],
        bed = config['anno_bed']
    shell:
        '''
        MisER -c {threads} \
        --strandSpecific \
        --setTag \
        --outBam {output} \
        {input} \
        {params.ref} \
        {params.bed} \
        {params.out_region} 2> {log}
        '''

# Replace "M" with "=" in CIGAR strings to generate PacBio compatible bam
rule modify_bam:
    input: "results/pbmm2/{id}.dedup.mapped.miser.corrected.bam"
    output: "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.bam"
    conda: "../envs/samtools_env.yaml"
    log: "results/logs/{id}_modify_bam.log"
    benchmark: "results/benchmarks/{id}_modify_bam.benchmark"
    threads: 16
    params:
        sam = "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sam"
    shell:
        '''
        samtools index -@ {threads} {input};
        samtools view -H {input} > {params.sam};
        samtools view {input} | awk 'BEGIN {{OFS="\\t"}} {{$6 = gensub("M", "=", "g", $6); print}}' >> {params.sam};
        sed 's/rc:i:0/rc:i:1/g' {params.sam} | samtools sort -@ {threads} -O BAM -o {output};
        samtools index -@ {threads} {output};
        rm {params.sam}
        '''

# Collapse reduntant transcripts into unique isoforms
# Alternative tool for this step would be [TAMA](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) or [Cupcake](https://github.com/Magdoll/cDNA_Cupcake) collapse script 
rule collapse:
    input: "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.bam" if config['is_pacbio'] else "results/pbmm2/{id}.dedup.mapped.miser.corrected.modified.sorted.ont.bam"
    output: "results/collapse/{id}.collapsed.gff"
    conda: "../envs/isoseq_env.yaml"
    benchmark: "results/benchmarks/{id}_collapse.benchmark"
    log: "results/logs/{id}_collapse.log"
    threads: 32
    shell:
        '''
        isoseq collapse -j {threads} \
        --log-file {log} \
        --log-level TRACE \
        {input} \
        {output}
        '''
