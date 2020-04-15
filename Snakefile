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
               mate=[1,2])

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

rule multiqc:
    output:
        "qc/multiqc_fastq.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc -n multiqc_fastq -o qc/ qc/fastq
        """
