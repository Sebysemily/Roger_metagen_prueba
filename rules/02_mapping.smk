import pandas as pd

configfile: "config/config.yml"

SAMPLES = pd.read_csv("config/samples.tsv", sep="\t")["sample"].tolist()

FILT_DIR     = config["paths"]["filtered_dir"]
HTML_DIR     = config["paths"]["html_dir"]
REF_DIR      = config["paths"]["ref_dir"]
MAPPING_DIR  = config["paths"]["mapping_dir"]
QC_DIR       = config["paths"]["qc_dir"]
DEHOSTED_DIR = config["paths"]["dehosted_dir"]

MAPQ_MIN = config["mapping"]["mapq_min"]
THREADS  = config["mapping"]["threads"]
EXTRA    = config["mapping"]["extra_options"]
HOSTS    = config["mapping"]["hosts"]   # ["mutuki","grandis","human"]

HOST_REF_FASTA = REF_DIR + "/hosts_combined.fna"
HOST_REF_INDEX = REF_DIR + "/hosts_combined.mmi"

# fija el tag => elimina wildcard {mapq} que rompía el DAG
MAPQ_TAG = f"mapq{MAPQ_MIN}"

# -------------------------------------------------------------------------
# 1) Concatenar genomas host en un solo FASTA
# -------------------------------------------------------------------------
rule build_host_reference:
    input:
        genomes=expand(REF_DIR + "/{host}.fna", host=HOSTS)
    output:
        ref=HOST_REF_FASTA
    conda:
        "../envs/02_mapping.yml"
    shell:
        r"""
        cat {input.genomes} > {output.ref}
        """

# -------------------------------------------------------------------------
# 2) Indexar el FASTA combinado
# -------------------------------------------------------------------------
rule index_host_reference:
    input:
        ref=HOST_REF_FASTA
    output:
        index=HOST_REF_INDEX
    conda:
        "../envs/02_mapping.yml"
    threads: THREADS
    shell:
        r"""
        minimap2 -x map-ont -d {output.index} {input.ref}
        """

# -------------------------------------------------------------------------
# 3) Mappear cada sample al host combinado
# -------------------------------------------------------------------------
rule map_sample_to_hostref:
    input:
        reads=FILT_DIR + "/{sample}_post_qc.fastq.gz",
        index=HOST_REF_INDEX
    output:
        bam=MAPPING_DIR + "/bam/{sample}.host.bam"
    conda:
        "../envs/02_mapping.yml"
    threads: THREADS
    shell:
        r"""
        mkdir -p {MAPPING_DIR}/bam
        minimap2 -t {threads} -ax map-ont {EXTRA} {input.index} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# -------------------------------------------------------------------------
# 4) Calcular estadisticas de mappeo (flagstat)
# -------------------------------------------------------------------------
rule flagstat:
    input:
        bam=MAPPING_DIR + "/bam/{sample}.host.bam"
    output:
        stats=MAPPING_DIR + "/flagstat/{sample}.host.flagstat"
    conda:
        "../envs/02_mapping.yml"
    shell:
        r"""
        mkdir -p {MAPPING_DIR}/flagstat
        samtools flagstat {input.bam} > {output.stats}
        """

# -------------------------------------------------------------------------
# 5) Generar CSV con % de mapeo
# -------------------------------------------------------------------------
rule mapping_stats:
    input:
        flagstats=expand(MAPPING_DIR + "/flagstat/{sample}.host.flagstat", sample=SAMPLES),
        nanostats=expand(QC_DIR + "/post/{sample}/NanoStats.txt", sample=SAMPLES)
    output:
        csv=MAPPING_DIR + "/stats_csv/mapping_stats.csv"
    script:
        "../code/mapping_stats.py"

# -------------------------------------------------------------------------
# 6) Extraer IDs mapeados al host con MAPQ >= MAPQ_MIN
# -------------------------------------------------------------------------
rule extract_mapped_ids:
    input:
        bam=MAPPING_DIR + "/bam/{sample}.host.bam"
    output:
        ids=MAPPING_DIR + f"/mapped_ids/{{sample}}.host_{MAPQ_TAG}.txt"
    conda:
        "../envs/02_mapping.yml"
    shell:
        r"""
        mkdir -p {MAPPING_DIR}/mapped_ids
        samtools view -F 4 -q {MAPQ_MIN} {input.bam} | cut -f1 | sort -u > {output.ids}
        """

# -------------------------------------------------------------------------
# 7) Quitar host (dehosted): filtra FASTQ original usando blacklist
# -------------------------------------------------------------------------
rule dehost_fastq:
    input:
        fastq=FILT_DIR + "/{sample}_post_qc.fastq.gz",
        blacklist=MAPPING_DIR + f"/mapped_ids/{{sample}}.host_{MAPQ_TAG}.txt"
    output:
        clean=DEHOSTED_DIR + "/{sample}_dehosted.fastq.gz"
    conda:
        "../envs/02_mapping.yml"
    shell:
        r"""
        mkdir -p {DEHOSTED_DIR}
        seqkit grep -v -f {input.blacklist} {input.fastq} | gzip > {output.clean}
        """

# -------------------------------------------------------------------------
# 8) Guardar REMOVED (lo que sí mapeó con MAPQ >= MAPQ_MIN)
# -------------------------------------------------------------------------
rule removed_fastq:
    input:
        fastq=FILT_DIR + "/{sample}_post_qc.fastq.gz",
        blacklist=MAPPING_DIR + f"/mapped_ids/{{sample}}.host_{MAPQ_TAG}.txt"
    output:
        removed=DEHOSTED_DIR + "/{sample}_removed.fastq.gz"
    conda:
        "../envs/02_mapping.yml"
    shell:
        r"""
        mkdir -p {DEHOSTED_DIR}
        seqkit grep -f {input.blacklist} {input.fastq} | gzip > {output.removed}
        """

# -------------------------------------------------------------------------
# 9) NanoPlot post-host
# -------------------------------------------------------------------------
rule Nanoplot_post_host:
    input:
        fastq=DEHOSTED_DIR + "/{sample}_dehosted.fastq.gz"
    output:
        post_host_stats=QC_DIR + "/post_host/{sample}/NanoStats.txt"
    conda:
        "../envs/01_qc.yml"
    threads: THREADS
    shell:
        r"""
        mkdir -p {QC_DIR}/post_host/{wildcards.sample}
        NanoPlot --fastq {input.fastq} -o {QC_DIR}/post_host/{wildcards.sample} \
          --plots dot --format png --dpi 200 -t {threads}
        """

# -------------------------------------------------------------------------
# 10) stats bases removidas
# -------------------------------------------------------------------------
rule bases_stats:
    input:
        nanostats_post=expand(QC_DIR + "/post/{sample}/NanoStats.txt", sample=SAMPLES),
        nanostats_post_host=expand(QC_DIR + "/post_host/{sample}/NanoStats.txt", sample=SAMPLES)
    output:
        csv=MAPPING_DIR + "/stats_csv/bases_stats.csv"
    conda:
        "../envs/02_mapping.yml"
    script:
        "../code/bases_stats.py"

# -------------------------------------------------------------------------
# 11) render HTML
# -------------------------------------------------------------------------
rule render_mapping_report:
    input:
        mapping_stats=MAPPING_DIR + "/stats_csv/mapping_stats.csv",
        bases_stats=MAPPING_DIR + "/stats_csv/bases_stats.csv"
    output:
        html=HTML_DIR + "/02_mapping.html"
    conda:
        "../envs/render.yml"
    shell:
        r"""
        mkdir -p {HTML_DIR}
        Rscript --vanilla -e '
        rmarkdown::render(
          "analysis/02_mapping.Rmd",
          params = list(
            mapping_stats = "{input.mapping_stats}",
            bases_stats   = "{input.bases_stats}"
          ),
          knit_root_dir = getwd()
        )
        '
        mv analysis/02_mapping.html {HTML_DIR}/02_mapping.html
        """
