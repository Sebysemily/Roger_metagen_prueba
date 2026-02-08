
import pandsas as pd

configfile: "config/config.yml"

SAMPLES = pd.read_csv("config/samples.tsv, sep="\t")["sample"].tolist()

RAW_DIR = config["paths"]["raw_dir] 
FILT_DIR = config["paths"]["filtered_dir]
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
	pre_ybl="{QC_DIR}/pre/{sample}/Yield_By_Length.html",
	pre_lvq="{QC_DIR}/pre/{sample}/LengthvsQualityScatterPlot_kde.html"
    conda: 
	"envs/01_qc.yml"
    threads: THREADS 
    shell: 
	"""
	mkdir -p {QC_DIR}/pre/{wildcards.sample}
	NanoPlot --fastq {input.fastq} -o {QC_DIR}/qc/pre/{wildcards.sample} --threads {threads} 
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

rule nanoplot_post:
    input:
        fastq=lambda wc: f"{FILT_DIR}/{wc.sample}_post_qc.fast.gz"
    output:
        post_stats="{QC_DIR}/post/{sample}/NanoStats.txt",
        post_ybl="{QC_DIR}/post/{sample}/Yield_By_Length.html",
        post_lvq="{QC_DIR}/post/{sample}/LengthvsQualityScatterPlot_kde.html"
    conda:
        "envs/01_qc.yml"
    threads: THREADS
    shell:
        """
        mkdir -p {QC_DIR}/post/{wildcards.sample}
        NanoPlot --fastq {input.fastq} -o {QC_DIR}/qc/post/{wildcards.sample} --threads {threads}
        """
# ---------------- REPORTE --------------------------------------------------------------------------
rule render_qc_report:
    input:
	pre=expand(f"{QC_DIR}/pre/{{sampe}}/NanoStats.txt}", sample=SAMPLES),
	post=expand(f"{QC_DIR}/post/{{sample}}/NanoStats.txt" sample=SAMPLES)
    output:
	hotml="results/results_docs/01_qc.html"
    shell:
	"""
	mkdir -p results/results_docs?01_qc.html
	"""
	
