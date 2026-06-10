SAMPLES = glob_wildcards('data_February_2026/{S}_R1_001.fastq.gz').S
READS = ['R1', 'R2']

rule all:
        input:
                expand("QC/pretrimming/{sample}_{read}_001_fastqc.zip", sample = SAMPLES, read = READS),
                expand("QC/pretrimming/{sample}_{read}_001_fastqc.html", sample = SAMPLES, read = READS)

rule fastqc_pre:
        input:
                "data_February_2026/{sample}_{read}_001.fastq.gz"
        output:
                "QC/pretrimming/{sample}_{read}_001_fastqc.zip",
                "QC/pretrimming/{sample}_{read}_001_fastqc.html"
        shell:
                "fastqc --outdir QC/pretrimming {input}"

