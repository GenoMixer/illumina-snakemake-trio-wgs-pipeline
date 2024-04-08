# libraries
import os
import glob
import numpy
import pandas as pd
from pathlib import Path

# config
#configfile: "config.yaml"
#WORK  = os.getcwd()
#head, tail = os.path.split(WORK)
#RUN = tail.split("_")[0]

# samples
samples=pd.read_csv("SampleSheet.csv").set_index("Sample_ID")
print(samples)

# variables
READ=["1", "2"]

# projet
project="240318_LH00212_0013_B222HMKLT1"

# functions
def get_fq1(wildcards):
        return sorted(glob.glob("raw_data/" + wildcards.sample + "*_R1_001.fastq.gz"))

def get_fq2(wildcards):
        return sorted(glob.glob("raw_data/" + wildcards.sample + "*_R2_001.fastq.gz"))

rule all:
    input:
         expand("results/{sample}/{sample}.sam", sample=samples.index),
         expand("results/{sample}/{sample}_sorted.bam", sample=samples.index),
         expand("results/{sample}/{sample}_sorted_marked.bam", sample=samples.index),
         expand("results/{sample}/{sample}_metrics.txt", sample=samples.index),
         expand("results/{sample}/{sample}_wgsmetrics.txt", sample=samples.index),
         expand("results/{sample}/{sample}_insertmetrics.txt", sample=samples.index),
         expand("results/{sample}/{sample}_inserthistogram.pdf", sample=samples.index),
         expand("results/{sample}/{sample}.grp", sample=samples.index),
         expand("results/{sample}/{sample}_improved.bam", sample=samples.index),
         expand("results/{sample}/{sample}.g.vcf.gz", sample=samples.index),
         expand("results/{project}/{project}.vcf.gz", project=project),
         expand("results/{project}/{project}.g.vcf.gz", project=project),
         expand("results/{sample}/{sample}.vcf.gz", sample=samples.index),
         expand("results/{sample}/{sample}.vcf.gz.tbi", sample=samples.index),
         expand("results/{sample}/{sample}.vcf.gz.md5", sample=samples.index),
         expand("results/{sample}/{sample}_varianteval.txt", sample=samples.index),
         expand("results/{sample}/{sample}_{read}_fastqc.html", sample=samples.index, read=READ),
         expand("results/{sample}/{sample}_{read}_fastqc.zip", sample=samples.index, read=READ),
         expand("results/QC/multiqc_report.html", sample=samples.index, read=READ)


rule BWA:
    input:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        fq1="raw_data/{sample}_1.fq.gz",
        fq2="raw_data/{sample}_2.fq.gz"
    output: 
        sam=temp("results/{sample}/{sample}.sam")
    params: 
        rg="@RG\\tID:{sample}\\tPL:Illumina\\tSM:{sample}\\tLB:WES"
    threads: 8
    log: 
        "logs/bwa/{sample}_bwa.log"
    resources:
        runtime="",
        mem_mb=""
    shell: 
        """
        bwa mem -R '{params.rg}' -t {threads} {input} > {output} 2> {log}
        """


rule SortSam:
    input:
        sam="results/{sample}/{sample}.sam"
    output:
        bam=temp("results/{sample}/{sample}_sorted.bam")
    params:
        sort_order="coordinate"
    log:
        "logs/picard/{sample}_sortsam.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        picard SortSam INPUT={input.sam} OUTPUT={output.bam} SORT_ORDER={params.sort_order} CREATE_INDEX=true CREATE_MD5_FILE=true 2> {log}
        """


rule MarkDuplicates:
    input: 
        bam="results/{sample}/{sample}_sorted.bam"
    output:
         bam=temp("results/{sample}/{sample}_sorted_marked.bam"),
         metrics=temp("results/{sample}/{sample}_metrics.txt")
    params: 
         extra="VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true CREATE_MD5_FILE=true"
    log: 
         "logs/picard/{sample}_markduplicates.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} {params.extra} 2> {log}
        """


rule BaseRecalibrator:
    input:
        bam="results/{sample}/{sample}_sorted_marked.bam",
    output:
        recal_table="results/{sample}/{sample}.grp",
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz",
        mills="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        indels="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
    log: 
        "logs/gatk/{sample}_baserecalibrator.log",
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk BaseRecalibrator --input {input.bam} --reference {params.ref} --known-sites {params.dbsnp} --known-sites {params.mills} --known-sites {params.indels} --output {output}  2> {log}
        """


rule ApplyBQSR:
    input:
        bam="results/{sample}/{sample}_sorted_marked.bam",
        recaltable="results/{sample}/{sample}.grp"
    output:
        "results/{sample}/{sample}_improved.bam"
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta"
    log:
        "logs/gatk/{sample}_applybqsr.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk ApplyBQSR --input {input.bam} --bqsr-recal-file {input.recaltable} --reference {params.ref} --output {output} 2> {log}
        """


rule HaplotypeCaller:
    input:
        bam="results/{sample}/{sample}_improved.bam",
    output:
        gvcf=temp("results/{sample}/{sample}.g.vcf.gz")
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta"
    log:
        "logs/gatk/{sample}_haplotypecaller.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk HaplotypeCaller -R {params.ref} -I {input.bam} -O {output.gvcf} -ERC GVCF 2> {log}
        """


rule CombineGVCF:
    input:
        gvcf=expand("results/{sample}/{sample}.g.vcf.gz", sample=samples.index),
    output:
        gvcf=temp(os.path.join("results/" + project + "/" + project + ".g.vcf.gz"))
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        gvcfs=lambda wildcards, input: list(map("--variant {}".format, input.gvcf))
    log:
        os.path.join("logs/gatk/" + project + "_combinegvcf.log")
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk CombineGVCFs {params.gvcfs} -R {params.ref} -O {output.gvcf}  2> {log}
        """


rule GenotypeGVCF:
    input:
        gvcf=os.path.join("results/" + project + "/" + project +  ".g.vcf.gz")
    output:
        vcf=os.path.join("results/" + project + "/" + project +  ".vcf.gz")
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz"
    log:
        os.path.join("logs/gatk/" + project + "_genotypegvcf.log")
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk GenotypeGVCFs -V {input.gvcf} -R {params.ref} -D {params.dbsnp} -O {output.vcf} 2> {log}
        """


rule GatherVCFs:
    input:
        vcf=os.path.join("results/" + project + "/" + project +  ".vcf.gz")
    output:
        vcf=protected("results/{sample}/{sample}.vcf.gz"),
        idx=protected("results/{sample}/{sample}.vcf.gz.tbi"),
        md5=protected("results/{sample}/{sample}.vcf.gz.md5")
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        extra=lambda wildcards, input: list(map("-sn ".format, wildcards.sample))
    log: 
        "logs/gatk/{sample}_gathervcf.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk SelectVariants -R {params.ref} -V {params.extra} -O {output} --create-output-variant-index --create-output-variant-md5 2> {log}
        """


rule FastQC:
    input:
        fastq="raw_data/{sample}_{read}.fq.gz"
    output: 
        temp("results/{sample}/{sample}_{read}_fastqc.zip"),
        temp("results/{sample}/{sample}_{read}_fastqc.html")
    log:
        "logs/fastqc/{sample}_{read}_fastqc.log"
    params:
        extra=""
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        fastqc {input} {output} 2> {log}
        """


rule CollectWgsMetrics:
    input:
        bam="results/{sample}/{sample}_improved.bam",
    output:
        metrics=temp("results/{sample}/{sample}_wgsmetrics.txt")
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        extra="",
        jvm_args=""
    log: 
        "logs/picard/{sample}_collect_wgsmetrics.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        picard CollectWgsMetrics INPUT={input.bam} OUTPUT={output} REFERENCE_SEQUENCE={params.ref}  2> {log}
        """


rule CollectInsertSizeMetrics:
    input:
        bam="results/{sample}/{sample}_improved.bam"
    output:
        metrics=temp("results/{sample}/{sample}_insertmetrics.txt"),
        histogram="results/{sample}/{sample}_inserthistogram.pdf"
    params:
        extra="M=0.5",
        jvm_args=""
    log: "logs/picard/{sample}_collect_insertsize_metrics.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        picard CollectInsertSizeMetrics INPUT={input.bam} OUTPUT={output.metrics} HISTOGRAM_FILE={output.histogram} 2> {log}
        """


rule VariantEval:
    input:
        vcf="results/{sample}/{sample}.vcf.gz",
    output: 
        eval="results/{sample}/{sample}_varianteval.txt"
    params:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz"
    log: 
        "logs/gatk/{sample}_varianteval.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk VariantEval -R {params.ref} -D {params.dbsnp} --eval {input.vcf} -O {output}
        """    


rule MultiQC:
    input:
        expand("results/{sample}/{sample}_{read}_fastqc.html", sample=samples.index, read=READ),
        expand("results/{sample}/{sample}_{read}_fastqc.zip", sample=samples.index, read=READ),
        expand("results/{sample}/{sample}_wgsmetrics.txt", sample=samples.index),
        expand("results/{sample}/{sample}_insertmetrics.txt", sample=samples.index),
        expand("results/{sample}/{sample}_metrics.txt", sample=samples.index),
        expand("results/{sample}/{sample}.grp", sample=samples.index),
        expand("results/{sample}/{sample}_varianteval.txt", sample=samples.index)
    output:
        "results/QC/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        multiqc -f -o unaligned {input} --no-data-dir
        """

