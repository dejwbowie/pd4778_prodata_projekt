configfile: "config.yaml"

SAMPLES = config["samples"]
READS = config["reads"]

RAW_DIR = config["raw_dir"]
RESULTS_DIR = config["results_dir"]
TRIMMED_DIR = f"{RESULTS_DIR}/trimmed"
FASTQC_DIR = f"{RESULTS_DIR}/fastqc"
MULTIQC_DIR = f"{RESULTS_DIR}/multiqc"
STAR_INDEX = "data/star_index"
ALIGNMENTS_DIR = f"{RESULTS_DIR}/alignments"
COVERAGE_DIR = f"{RESULTS_DIR}/coverage"
VARIANTS_DIR = f"{RESULTS_DIR}/variants"
REF_GENOME = config["ref_genome"]

rule all:
    input:
        expand(f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=READS),
        expand(f"{TRIMMED_DIR}/{{sample}}_1_trimmed.fastq", sample=SAMPLES),
        expand(f"{TRIMMED_DIR}/{{sample}}_2_trimmed.fastq", sample=SAMPLES),
        expand(f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(f"{COVERAGE_DIR}/{{sample}}.bedgraph", sample=SAMPLES),
        expand(f"{VARIANTS_DIR}/{{sample}}.vcf", sample=SAMPLES),
        f"{MULTIQC_DIR}/multiqc_report.html"

rule prepare_dirs:
    output:
        directory(RESULTS_DIR)
    run:
        import os
        os.makedirs(FASTQC_DIR, exist_ok=True)
        os.makedirs(MULTIQC_DIR, exist_ok=True)
        os.makedirs(TRIMMED_DIR, exist_ok=True)
        os.makedirs(ALIGNMENTS_DIR, exist_ok=True)
        os.makedirs(COVERAGE_DIR, exist_ok=True)
        os.makedirs(VARIANTS_DIR, exist_ok=True)
        os.makedirs(STAR_INDEX, exist_ok=True)

rule fastqc:
    input:
        lambda wc: f"{RAW_DIR}/{wc.sample}_{wc.read}.fastq"
    output:
        html = temp(f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.html"),
        zip = f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.zip"
    log:
        f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.log"
    container:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: 2
    shell:
        "fastqc {input} --outdir {FASTQC_DIR} --threads {threads} > {log} 2>&1"

rule multiqc:
    input:
        expand(f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=READS)
    output:
        f"{MULTIQC_DIR}/multiqc_report.html"
    container:
        "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    shell:
        "multiqc {FASTQC_DIR} -o {MULTIQC_DIR}"

rule trimmomatic:
    input:
        r1 = lambda wc: f"{RAW_DIR}/{wc.sample}_1.fastq",
        r2 = lambda wc: f"{RAW_DIR}/{wc.sample}_2.fastq"
    output:
        trimmed_r1 = f"{TRIMMED_DIR}/{{sample}}_1_trimmed.fastq",
        trimmed_r2 = f"{TRIMMED_DIR}/{{sample}}_2_trimmed.fastq",
        unpaired_r1 = temp(f"{TRIMMED_DIR}/{{sample}}_1_unpaired.fastq"),
        unpaired_r2 = temp(f"{TRIMMED_DIR}/{{sample}}_2_unpaired.fastq")
    container:
        "docker://quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    threads: 4
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.trimmed_r1} {output.unpaired_r1} \
            {output.trimmed_r2} {output.unpaired_r2} \
            ILLUMINACLIP:/usr/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule star_index:
    input:
        ref = REF_GENOME
    output:
        touch(f"{STAR_INDEX}/SA")
    container:
        "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_1"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {STAR_INDEX} \
            --genomeFastaFiles {input.ref}
        """

rule star_align:
    input:
        r1 = lambda wc: f"{TRIMMED_DIR}/{wc.sample}_1_trimmed.fastq",
        r2 = lambda wc: f"{TRIMMED_DIR}/{wc.sample}_2_trimmed.fastq",
        index = f"{STAR_INDEX}/SA"
    output:
        bam = f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam"
    container:
        "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_1"
    threads: 8
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {STAR_INDEX} \
            --readFilesIn {input.r1} {input.r2} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {ALIGNMENTS_DIR}/{wildcards.sample}_
        """

rule index_bam:
    input:
        f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam"
    output:
        f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    container:
        "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    shell:
        "samtools index {input}"

rule coverage:
    input:
        bam = f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    output:
        f"{COVERAGE_DIR}/{{sample}}.bedgraph"
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    shell:
        "bedtools genomecov -ibam {input.bam} -bg > {output}"

rule variant_calling:
    input:
        bam = f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{ALIGNMENTS_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
        ref = REF_GENOME
    output:
        f"{VARIANTS_DIR}/{{sample}}.vcf"
    container:
        "docker://quay.io/biocontainers/freebayes:1.3.5--h2d02072_3"
    shell:
        "freebayes -f {input.ref} {input.bam} > {output}"
