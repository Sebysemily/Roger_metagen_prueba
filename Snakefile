import pandas as pd

configfile: "config/config.yml"

SAMPLES = pd.read_csv("config/samples.tsv", sep="\t")["sample"].tolist()

HTML_DIR     = config["paths"]["html_dir"]
DEHOSTED_DIR = config["paths"]["dehosted_dir"]
BLAST_DIR    = config["paths"]["blast_dir"]

include: "rules/01_qc.smk"
include: "rules/02_mapping.smk"
include: "rules/03_blast.smk"

rule all:
    input:
        # Reporte QC
        HTML_DIR + "/01_qc.pdf",

        # Reporte mapping
        HTML_DIR + "/02_mapping.html",

        # Outputs que 03_blast necesita (generados en 02_mapping)
        expand(DEHOSTED_DIR + "/{sample}_removed.fastq.gz",  sample=SAMPLES),
        expand(DEHOSTED_DIR + "/{sample}_dehosted.fastq.gz", sample=SAMPLES),

        # Outputs finales de BLAST
        expand(BLAST_DIR + "/local_hosts/{sample}.removed.vs_hosts.tsv",   sample=SAMPLES),
        expand(BLAST_DIR + "/remote_nt/{sample}.dehosted.vs_remote.tsv", sample=SAMPLES),
        expand(BLAST_DIR + "/remote_nt/{sample}.dehosted.vs_remote.summary.csv", sample=SAMPLES)
