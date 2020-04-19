configfile: "config.yaml"

import pandas as pd
import os

# set display options
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)


metadata = pd.read_csv(config["metadataFile"], dtype=str, sep='\t',
                       usecols=['Source Name', 'Comment[FASTQ_URI]']).rename(columns={
                           'Source Name': 'sample_id',
                           'Comment[FASTQ_URI]': 'url'})
metadata['mate'] = "R" + metadata.url.str.extract(r'(1|2)\.fastq')
metadata['cell_type'] = metadata.sample_id.str.extract(r'(HSC|MPP[1234])')
metadata['replicate'] = metadata.sample_id.str.extract(r'(.)$')

metadata = metadata.pivot_table(index=["sample_id", "cell_type", "replicate"],
                                columns="mate", values="url", aggfunc=lambda x: x[0])
metadata = metadata.reset_index().set_index('sample_id')
#print(metadata)


rule all:
    input:
        expand("qc/fastq/{sample_id}_{mate}_fastqc.html",
               sample_id=metadata.index.tolist(),
               mate=[1,2]),
        "data/se/hspcs_bulk.txs.rds",
        "data/se/hspcs_bulk.genes.rds"

rule download_fastq:
    output:
        r1="data/fastq/{sample_id}_1.fastq.gz",
        r2="data/fastq/{sample_id}_2.fastq.gz"
    params:
        url1=lambda wcs: metadata.R1[wcs.sample_id],
        url2=lambda wcs: metadata.R2[wcs.sample_id]
    shell:
        """
        wget -O {output.r1} {params.url1}
        wget -O {output.r2} {params.url2}
        """

rule fastqc:
    input:
        lambda wcs: "data/fastq/%s_%s.fastq.gz" % (wcs.sample_id, wcs.mate)
    output:
           "qc/fastq/{sample_id}_{mate}_fastqc.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        outDir=$(dirname {output})
        fastqc -o $outDir {input}
        """

rule multiqc_fastq:
    output:
        "qc/multiqc_fastq.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc -n multiqc_fastq -o qc/ qc/fastq
        """

rule hisat2_r2:
    input:
        "data/fastq/{sample_id}_2.fastq.gz"
    output:
        bam="data/bam/hisat2_r2/{sample_id}.bam",
        bai="data/bam/hisat2_r2/{sample_id}.bam.bai"
    params:
        tmp=config['tmp_dir'],
        hisat2=config['hisat2'],
        idx=config['hisat2Index'],
        samtools=config['samtools'],
        sam=config['tmp_dir'] + "/{sample_id}.sam"
    threads: 16
    resources:
        mem=2
    shell:
        """
        {params.hisat2} -p {threads} -x {params.idx} -U {input} -S {params.sam}
        {params.samtools} sort -@ {threads} -m 2G -T {params.tmp}/ -o {output.bam} {params.sam}
        {params.samtools} index -@ {threads} {output.bam}
        rm -f {params.sam}
        """

rule kallisto_quant_r2:
    input:
        fastq="data/fastq/{sample_id}_2.fastq.gz",
        kdx=config['utromeKDX'],
        gtf=config['utromeGTF'],
        chroms=config['chromSizes']
    output:
        bam="data/bam/utrome_r2/{sample_id}.bam",
        bai="data/bam/utrome_r2/{sample_id}.bam.bai",
        tsv="data/kallisto/utrome_r2/{sample_id}/abundance.tsv"
    conda:
        "envs/kallisto.yaml"
    threads: 8
    resources:
        mem=1
    shell:
        """
        outDir=$(dirname {output.tsv})
        rm -rf $outDir
        kallisto quant -t {threads} -i {input.kdx} \\
        --genomebam -g {input.gtf} -c {input.chroms} \\
        --bias --single --single-overhang --fr-stranded -l1 -s1 \\
        -o $outDir {input.fastq}
        mv $outDir/pseudoalignments.bam {output.bam}
        mv $outDir/pseudoalignments.bam.bai {output.bai}
        """

rule multiqc_kallisto:
    output:
        "qc/multiqc_kallisto.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc -n multiqc_kallisto -o qc/ logs/lsf/kallisto_quant*
        """

rule kallisto_to_se:
    output:
        "data/se/hspcs_bulk.txs.raw.rds"
    input:
        scriptFile="scripts/kallisto_to_se.R",
        tx2gene=config['utromeTx2Gene'],
        tsvs=expand("data/kallisto/utrome_r2/{cell_type}{batch}/abundance.h5",
                    cell_type=["HSC", "MPP1", "MPP2", "MPP3", "MPP4"],
                    batch=["a", "b", "c", "d"])
    params:
        inputDir="data/kallisto/utrome_r2"
    resources:
        mem=1,
        walltime=1
    conda:
        "envs/r36-tximport.yaml"
    shell:
        """
        {input.scriptFile} {params.inputDir} {input.tx2gene} {output}
        """

rule collapse_overlaps:
    input:
        scriptFile="scripts/collapse_se_txs.R",
        se="data/se/hspcs_bulk.txs.raw.rds",
        txMapFile=config['utromeMerge'],
        annotsFile=config['atlasUTRAnnots']
    output:
        "data/se/hspcs_bulk.txs.rds"
    resources:
        mem=2,
        walltime=1
    conda:
        "envs/r36-tximport.yaml"
    shell:
        """
        {input.scriptFile} {input.se} {input.txMapFile} {input.annotsFile} {output}
        """

rule se_txs_to_genes:
    input:
        scriptFile="scripts/se_txs_to_genes.R",
        se="data/se/hspcs_bulk.txs.rds",
        annots=config['atlasGeneAnnots']
    output:
        "data/se/hspcs_bulk.genes.rds"
    resources:
        mem=2,
        walltime=1
    conda:
        "envs/r36-tximport.yaml"
    shell:
        """
        {input.scriptFile} {input.se} {input.annots} {output}
        """
