
import pandsas as pd

configfile: "config/config.yml"

SAMPLES = pd.read_csv("config/samples.tsv", sep="\t")["sample"].tolist()

RAW_DIR = config["paths"]["raw_dir"] 
FILT_DIR = config["paths"]["filtered_dir"]
QC_DIR = config["paths"]["qc_dir"]

MIN_LEN = config["qc"]["min_length"]
MIN_Q = config["qc"]["min_mean_q"]
THREADS = config["qc"]["threads"]

rule all: 
    input:
        "results/results_docs/01_qc.html"
# -------------- QC PRE ------------------------------------------------------------------------------
rule nanoplot_pre:
    input:
        fastq=lambda wc: f"{RAW_DIR}/{wc.sample}.fast.gz"
    output:
        pre_stats="{QC_DIR}/pre/{sample}/NanoStats.txt",
        pre_ybl="{QC_DIR}/pre/{sample}/Yield_By_Length.png",
        pre_lvq="{QC_DIR}/pre/{sample}/LengthvsQualityScatterPlot_kde.png"
    conda: 
        "envs/01_qc.yml"
        threads: THREADS 
    shell: 
        """
        mkdir -p {QC_DIR}/pre/{wildcards.sample}
        NanoPlot --fastq {input.fastq} -o {QC_DIR}/qc/pre/{wildcards.sample} --plots dot --format png --dpi 200 -t {threads} 
        """
# ------------ FILTRADO --------------------------------------------------------------------------
rule fitlong:
    input:
        fastq=lambda wc: f"{RAW_DIR}/{wc.sample}.fast.gz"
    output:
        fastq_filt=lambda wc: f"{FILT_DIR}/{wc.sample}_post_qc.fastq.gz"	
    shell:
        """
        mkdir -p {FILT_DIR}
        filtlong --min_length {MIN_LEN} --min_mean_q {MIN_Q} | gzip > {output.fastq_filt}
        """
# ------------- QC POST ----------------------------------------------------------------
rule nanoplot_post:
    input:
        fastq=lambda wc: f"{FILT_DIR}/{wc.sample}_post_qc.fast.gz"
    output:
        post_stats="{QC_DIR}/post/{sample}/NanoStats.txt",
        post_ybl="{QC_DIR}/post/{sample}/Yield_By_Length.png",
        post_lvq="{QC_DIR}/post/{sample}/LengthvsQualityScatterPlot_kde.png"
    conda:
        "envs/01_qc.yml"
    threads: THREADS
    shell:
        """
        mkdir -p {QC_DIR}/post/{wildcards.sample}
        NanoPlot --fastq {input.fastq} -o {QC_DIR}/qc/post/{wildcards.sample} --plots dot --format png --dpi 200 -t {threads}
        """
# ---------------- REPORTE --------------------------------------------------------------------------
rule render_qc_report:
    input:
        pre_stats=expand(f"{QC_DIR}/pre/{{sample}}/NanoStats.txt", sample=SAMPLES),
        post_stats=expand(f"{QC_DIR}/post/{{sample}}/NanoStats.txt", sample=SAMPLES),
        pre_ybl=expand(f"{QC_DIR}/pre/{{sample}}/Yield_By_Length.png", sample=SAMPLES),
        post_ybl=expand(f"{QC_DIR}/post/{{sample}}/Yield_By_Length.png", sample=SAMPLES),   
        pre_lvq=expand(f"{QC_DIR}/pre/{{sample}}/LengthvsQualityScatterPlot_kde.png", sample=SAMPLES),
        post_lvq=expand(f"{QC_DIR}/post/{{sample}}/LengthvsQualityScatterPlot_kde.png", sample=SAMPLES) 
    output:
        html="results/results_docs/01_qc.html"
    shell:
        """
        mkdir -p results/results_docs
        Rscript --vanilla -e \
        'rmarkdown::render(
        "analysis/01_qc.Rmd",
        output_file = "{output.html}",
        params = list(
        pre_stats = "{input.pre_stats}",
        post_stats = "{input.post_stats}",
        pre_ybl = "{input.pre_ybl}",
        post_ybl = "{input.post_ybl}",
        pre_lvq = "{input.pre_lvq}",
        post_lvq = "{input.post_lvq}"
         )
        )'
        """ 

