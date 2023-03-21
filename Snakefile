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
reference="/mnt/ngs-resources/references/bwa/human_hs38DH/hs38DH.fa"

# samples
samples=pd.read_csv("SampleSheet.csv").set_index("Sample_ID")
print(samples)

# variables
READ=["1", "2"]

# functions
def get_fq1(wildcards):
        return sorted(glob.glob("raw_data/" + wildcards.sample + "*_R1_001.fastq.gz"))

def get_fq2(wildcards):
        return sorted(glob.glob("raw_data/" + wildcards.sample + "*_R2_001.fastq.gz"))

rule all:
    input:
         expand("raw_data/{sample}_1.fq.gz", sample=samples.index),
         expand("raw_data/{sample}_2.fq.gz", sample=samples.index),
         expand("raw_data/{sample}_qc.txt", sample=samples.index),
         expand("alignment/{sample}.sam", sample=samples.index),
         expand("alignment/{sample}_sorted.bam", sample=samples.index),
         expand("alignment/{sample}_sorted_marked.bam", sample=samples.index),
         expand("alignment/{sample}_metrics.txt", sample=samples.index),
         expand("qc/multiqc_report.html", sample=samples.index, read=READ)


rule cutadapt:
    input:
        fq1=get_fq1,
        fq2=get_fq2
    output:
        op1="raw_data/{sample}_1.fq.gz",
        op2="raw_data/{sample}_2.fq.gz",
        qc="raw_data/{sample}_qc.txt"
    params:
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        extra="-Z --overlap 8 --error-rate 0.1 --no-indels -m 30 -q 20,20 --trim-n"
    threads: 2
    log: "logs/cutadapt/{sample}_trim.log"
    shell: 
        """
        cutadapt --cores {threads} {params.adapters} {params.extra} -o {output.op1} -p {output.op2} {input} > {output.qc} 2> {log}
        """

rule bwa:
    input:
        ref="/mnt/ngs-resources/references/bwa/human_hs38DH/hs38DH.fa",
        fq1="raw_data/{sample}_1.fq.gz",
        fq2="raw_data/{sample}_2.fq.gz"
    output: "alignment/{sample}.sam"
    params: rg="@RG\\tID:{sample}\\tPL:Illumina\\tSM:{sample}\\tLB:WES"
    threads: 1
    log: "logs/bwa_mem/{sample}_aln.log"
    shell: 
        """
        bwa mem -R '{params.rg}' -t {threads} {input} > {output} 2> {log}
        """

rule sortsam:
    input: "alignment/{sample}.sam"
    output: "alignment/{sample}_sorted.bam"
    params: sort_order="coordinate"
    log: "logs/picard/{sample}_sortsam.log"
    shell:
        """
        picard SortSam INPUT={input} OUTPUT={output} SORT_ORDER={params.sort_order} CREATE_INDEX=true CREATE_MD5_FILE=true 2> {log}
        """ 	

rule markduplicates:
    input: "alignment/{sample}_sorted.bam"
    output: 
         bam="alignment/{sample}_sorted_marked.bam",
         metrics="alignment/{sample}_metrics.txt"
    params: extra="VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true CREATE_MD5_FILE=true"
    log: "logs/picard/{sample}_markduplicates.log"
    shell:
        """
        picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} {params.extra} 2> {log}
        """

rule fastqc:
    input: expand("raw_data/{sample}_{read}.fq.gz", sample=samples.index, read=READ)
    output: 
        temp("raw_data/{sample}_{read}_fastqc.zip"),
        temp("raw_data/{sample}_{read}_fastqc.html")
    shell: "fastqc {input} {output}"

rule multiqc:
     input: 
        expand("raw_data/{sample}_{read}_fastqc.html", sample=samples.index, read=READ),
        expand("raw_data/{sample}_{read}_fastqc.zip", sample=samples.index, read=READ),
        expand("raw_data/{sample}_{read}_qc.txt", sample=samples.index, read=READ),
        expand("alignment/{sample}_metrics.txt", sample=samples.index, read=READ)
     output: "qc/multiqc_report.html"
     shell: "multiqc -f -o unaligned {input} --no-data-dir"



