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

# references
reference="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta"

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
         expand("alignment/{sample}.sam", sample=samples.index),
         expand("alignment/{sample}_sorted.bam", sample=samples.index),
         expand("alignment/{sample}_sorted_marked.bam", sample=samples.index),
         expand("alignment/{sample}_metrics.txt", sample=samples.index),
         expand("alignment/{sample}_wgsmetrics.txt", sample=samples.index),
         expand("alignment/{sample}_insertmetrics.txt", sample=samples.index),
         expand("alignment/{sample}_inserthistogram.pdf", sample=samples.index),
         expand("alignment/{sample}.grp", sample=samples.index),
         expand("alignment/{sample}_improved.bam", sample=samples.index),
         expand("alignment/{sample}.vcf.gz", sample=samples.index),
#         expand("alignment/{project}.vcf.gz", project=project),
         expand("alignment/{sample}_varianteval.txt", sample=samples.index),
         expand("qc/multiqc_report.html", sample=samples.index, read=READ)

rule bwa:
    input:
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        fq1="raw_data/{sample}_1.fq.gz",
        fq2="raw_data/{sample}_2.fq.gz",
    output: 
        sam=temp("alignment/{sample}.sam")
    params: 
        rg="@RG\\tID:{sample}\\tPL:Illumina\\tSM:{sample}\\tLB:WES"
    threads: 1
    log: 
        "logs/bwa/{sample}_aln.log"
    shell: 
        """
        bwa mem -R '{params.rg}' -t {threads} {input} > {output} 2> {log}
        """

rule sortsam:
    input:
        "alignment/{sample}.sam"
    output:
        temp("alignment/{sample}_sorted.bam")
    params:
        sort_order="coordinate"
    log:
        "logs/picard/sortsam/{sample}_sortsam.log"
    shell:
        """
        picard SortSam INPUT={input} OUTPUT={output} SORT_ORDER={params.sort_order} CREATE_INDEX=true CREATE_MD5_FILE=true 2> {log}
        """

rule markduplicates:
    input: "alignment/{sample}_sorted.bam"
    output:
         bam=temp("alignment/{sample}_sorted_marked.bam"),
         metrics=temp("alignment/{sample}_metrics.txt"),
    params: extra="VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true CREATE_MD5_FILE=true"
    log: "logs/picard/markduplicated/{sample}_markduplicates.log"
    shell:
        """
        picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} {params.extra} 2> {log}
        """

rule baserecalibrator:
    input:
        bam="alignment/{sample}_sorted_marked.bam",
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz",
        mills="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        indels="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
    output:
        recal_table="alignment/{sample}.grp",
    log: 
        "logs/gatk/baserecalibrator/{sample}.log",
    shell:
        """
        gatk BaseRecalibrator --input {input.bam} --reference {input.ref} --known-sites {input.dbsnp} --known-sites {input.mills} --known-sites {input.indels} --output {output}  2> {log}
        """

rule applybqsr:
    input:
        bam="alignment/{sample}_sorted_marked.bam",
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        recaltable="alignment/{sample}.grp"
    output:
        "alignment/{sample}_improved.bam"
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    shell:
        """
        gatk ApplyBQSR --input {input.bam} --bqsr-recal-file {input.recaltable} --reference {input.ref} --output {output} 2> {log}
        """

rule haplotypecaller:
    input:
        bam="alignment/{sample}_improved.bam",
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta"
    output:
        gvcf=temp("alignment/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.gvcf} -ERC GVCF 2> {log}
        """

rule genotypegvcf:
    input:
        gvcf="alignment/{sample}.g.vcf.gz",
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz"
    output:
        vcf="alignment/{sample}.vcf.gz"
    log:
        "logs/gatk/genotypegvcf/{sample}.log"
    shell:
        """
        gatk GenotypeGVCFs -R {input.ref} -D {input.dbsnp} -V {input.gvcf} -O {output.vcf} 2> {log}
        """

#rule combinegvcfs:
#    input:
#        gvcfs=expand("alignment/{sample}.g.vcf.gz", sample=samples.index),
#    output:
#        vcf="alignment/{project}.vcf.gz"
#    params:
#        project="{project}",
#        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
#        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz",
#        extra=lambda wildcards, input: ' --variant '.join(input.gvcfs)
#    log:
#        "logs/gatk/combinegvcf/{project}.log"
#    shell:
#        """
#        gatk CombineGVCFs -R {params.ref} -D {params.dbsnp} {params.extra} -O {output.vcf}  2> {log}
#        """

rule fastqc:
    input:
        expand("raw_data/{sample}_{read}.fq.gz", sample=samples.index, read=READ)
    output: 
        temp("raw_data/{sample}_{read}_fastqc.zip"),
        temp("raw_data/{sample}_{read}_fastqc.html")
    shell:
        """
        fastqc {input} {output} 2> {log}
        """

rule collect_wgsmetrics:
    input:
        bam="alignment/{sample}_improved.bam",
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta"
    output:
        metrics=temp("alignment/{sample}_wgsmetrics.txt")
    params:
        extra="",
        jvm_args=""
    log: "logs/picard/{sample}_collect_wgsmetrics.log"
    resources:
        runtime="",
        mem_mb="",
    shell:
        """
        picard CollectWgsMetrics INPUT={input.bam} OUTPUT={output} REFERENCE_SEQUENCE={input.ref}  2> {log}
        """


rule collect_insertsize_metrics:
    input:
        bam="alignment/{sample}_improved.bam"
    output:
        metrics=temp("alignment/{sample}_insertmetrics.txt"),
        histogram="alignment/{sample}_inserthistogram.pdf"
    params:
        extra="M=0.5",
        jvm_args=""
    log: "logs/picard/{sample}_collect_insertsize_metrics.log"
    resources:
        runtime="",
        mem_mb="",
    shell:
        """
        picard CollectInsertSizeMetrics INPUT={input.bam} OUTPUT={output.metrics} HISTOGRAM_FILE={output.histogram} 2> {log}
        """

rule variant_eval:
    input:
        vcf="alignment/{sample}.vcf.gz",
        ref="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/Homo_sapiens_assembly38.fasta",
        dbsnp="/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/resources/genomes/human/gatk/hg38/dbsnp_138.hg38.vcf.gz"
    output: eval="alignment/{sample}_varianteval.txt"
    log: "logs/{sample}/variant_eval.log"
    resources:
        runtime="",
        mem_mb=""
    shell:
        """
        gatk VariantEval -R {input.ref} -D {input.dbsnp} --eval {input.vcf} -O {output}
        """    

rule multiqc:
    input:
#        expand("raw_data/{sample}_{read}_fastqc.html", sample=samples.index, read=READ),
#        expand("raw_data/{sample}_{read}_fastqc.zip", sample=samples.index, read=READ),
        expand("alignment/{sample}_wgsmetrics.txt", sample=samples.index),
        expand("alignment/{sample}_insertmetrics.txt", sample=samples.index),
        expand("alignment/{sample}_metrics.txt", sample=samples.index),
        expand("alignment/{sample}.grp", sample=samples.index),
        expand("alignment/{sample}_varianteval.txt", sample=samples.index)
    output:
        "qc/multiqc_report.html"
    shell:
        """
        multiqc -f -o unaligned {input} --no-data-dir
        """

