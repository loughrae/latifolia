SAMPLES = glob_wildcards('data_February_2026/{S}_R1_001.fastq.gz').S
READS = ['R1', 'R2']

DNAS = glob_wildcards('QC/posttrimming/VB{dnaS}_R1_001.trimmed_fastp.fastq.gz').dnaS

rule all:
        input:
                expand("QC/pretrimming/{sample}_{read}_001_fastqc.zip", sample = SAMPLES, read = READS),
                expand("QC/pretrimming/{sample}_{read}_001_fastqc.html", sample = SAMPLES, read = READS),
                expand("QC/posttrimming/{sample}_{read}_001.trimmed_fastp.fastq.gz", sample = SAMPLES, read = READS),
                expand("QC/posttrimming/{sample}_001.fastpqc.html", sample = SAMPLES),
                expand("QC/posttrimming/{sample}_001.fastpqc.json", sample = SAMPLES),
                expand("QC/posttrimming/{sample}_{read}_001.trimmed_fastp_fastqc.zip", sample = SAMPLES, read = READS),
                expand("QC/posttrimming/{sample}_{read}_001.trimmed_fastp_fastqc.html", sample = SAMPLES, read = READS),
                expand("alignments/VB{dnas}.bam", dnas = DNAS),
                expand("alignments/VB{dnas}_markdup.bam", dnas = DNAS),
                expand("alignments/VB{dnas}_markdup.bam.bai", dnas = DNAS),
                expand("alignments/VB{dnas}.vcf.gz", dnas = DNAS),
                expand("alignments/VB{dnas}_markdup_RG.bam", dnas = DNAS),
                expand("alignments/VB{dnas}_markdup_RG.bam.bai", dnas = DNAS),
                expand("alignments/VB{dnas}_gatk.vcf.gz", dnas = DNAS)

rule fastqc_pre:
        input:
                "data_February_2026/{sample}_{read}_001.fastq.gz"
        output:
                "QC/pretrimming/{sample}_{read}_001_fastqc.zip",
                "QC/pretrimming/{sample}_{read}_001_fastqc.html"
        shell:
                "fastqc --outdir QC/pretrimming {input}"

rule fastp:
        input:
                in1="data_February_2026/{sample}_R1_001.fastq.gz",
                in2="data_February_2026/{sample}_R2_001.fastq.gz"
        output:
                out1="QC/posttrimming/{sample}_R1_001.trimmed_fastp.fastq.gz",
                out2="QC/posttrimming/{sample}_R2_001.trimmed_fastp.fastq.gz",
                out3="QC/posttrimming/{sample}_001.fastpqc.html",
                out4="QC/posttrimming/{sample}_001.fastpqc.json"
        shell:
                """
                fastp -i {input.in1} -I {input.in2} -o {output.out1} -O {output.out2} --trim_front1 15 \
                -h {output.out3} -j {output.out4}
                """

rule fastqc_post:
        input:
                "QC/posttrimming/{sample}_{read}_001.trimmed_fastp.fastq.gz"
        output:
                "QC/posttrimming/{sample}_{read}_001.trimmed_fastp_fastqc.zip",
                "QC/posttrimming/{sample}_{read}_001.trimmed_fastp_fastqc.html"
        shell:
                "fastqc --outdir QC/posttrimming {input}"


#I have already created the genome index: bwa-mem2 index Slatifolia.v5.genome.fna 

rule bwa_mem2:
        threads: 14
        input:
                in1="QC/posttrimming/VB{dnas}_R1_001.trimmed_fastp.fastq.gz",
                in2="QC/posttrimming/VB{dnas}_R2_001.trimmed_fastp.fastq.gz",
                ref="Silene_latifolia_genome_v5/Slatifolia.v5.genome.fna"
        output:
                "alignments/VB{dnas}.bam"
        shell:
                "bwa-mem2 mem -t 12 {input.ref} {input.in1} {input.in2} | samtools sort -@4 -o {output}"

rule mark_dups:
        input:
                "alignments/VB{dnas}.bam"
        output:
                "alignments/VB{dnas}_markdup.bam"
        shell:
                "samtools collate -Ou {input} | samtools fixmate -m -u - - | samtools sort -u - -o - | samtools markdup - {output}"

rule make_dna_index:
        input:
                "alignments/VB{dnas}_markdup.bam"
        output:
                "alignments/VB{dnas}_markdup.bam.bai"
        shell:
                "samtools index {input}"

rule make_vcfs:
        input:
               in1= "alignments/VB{dnas}_markdup.bam",
               ploidy="alignments/VB{dnas}.ploidy",
               ref="Silene_latifolia_genome_v5/Slatifolia.v5.genome.fna"
        output:
               "alignments/VB{dnas}.vcf.gz"
        shell:
               "bcftools mpileup -f {input.ref} -Ou {input.in1} | bcftools call --ploidy-file {input.ploidy} -mv -Oz -o {output}"

PLOIDIES = {
    "1_S71": 2,
    "2_S72": 2,
    "3_S73": 4,
}

rule add_read_groups:
    input:
        "alignments/VB{dnas}_markdup.bam"
    output:
        "alignments/VB{dnas}_markdup_RG.bam"
    conda:
        "gatk46"
    shell:
         "gatk AddOrReplaceReadGroups -I {input}  -O {output}  --RGID VB{wildcards.dnas} --RGLB lib1  --RGPL ILLUMINA --RGPU unit1  --RGSM VB{wildcards.dnas}"
       
rule make_dna_index_RG:
        input:
                "alignments/VB{dnas}_markdup_RG.bam"
        output:
                "alignments/VB{dnas}_markdup_RG.bam.bai"
        shell:
                "samtools index {input}"

rule make_vcfs_gatk:
        input:
               in1="alignments/VB{dnas}_markdup_RG.bam",
               ploidyregions="alignments/VB{dnas}.bed",
               ref="Silene_latifolia_genome_v5/Slatifolia.v5.genome.fna"
        params:
               ploidy=lambda wc: PLOIDIES[wc.dnas]
        conda:
               "gatk46"
        output:
               "alignments/VB{dnas}_gatk.vcf.gz"
        shell:
               "gatk HaplotypeCaller -R {input.ref} -I {input.in1} -ploidy {params.ploidy}  --ploidy-regions {input.ploidyregions} -O {output}"
        
