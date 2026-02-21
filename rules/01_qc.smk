
import pandas as pd
import glob
import os

configfile: "config/config.yml"

samples_info = pd.read_csv("config/samples.tsv", sep="\t")

SAMPLES = samples_info["sample"].tolist()
BARCODES = samples_info["barcode"].tolist()
sample_to_barcode = dict(zip(samples_info["sample"], samples_info["barcode"]))

RAW_DIR = config["paths"]["raw_dir"] 
FILT_DIR = config["paths"]["filtered_dir"]
QC_DIR = config["paths"]["qc_dir"]
HTML_DIR = config["paths"]["html_dir"]

MIN_LEN = config["qc"]["min_length"]
MIN_Q = config["qc"]["min_mean_q"]
THREADS = config["qc"]["threads"]
#--------------- CONC SAMPLES-----------------------------------------------
rule concatenate_sample:
    input:
        lambda wildcards: sorted(glob.glob(os.path.join(RAW_DIR, sample_to_barcode[wildcards.sample], "*.fastq.gz")))
    output:
        combined = RAW_DIR + "/combined/{sample}.fastq.gz"
    run:
        import gzip
        os.makedirs(os.path.dirname(output.combined), exist_ok=True)
        with gzip.open(output.combined, 'wb') as out:
            for f in input:
                with gzip.open(f, 'rb') as inf:
                    out.write(inf.read())
# -------------- QC PRE ------------------------------------------------------------------------------
rule nanoplot_pre:
    input:
        fastq=RAW_DIR + "/combined/{sample}.fastq.gz"
    output:
        pre_stats=QC_DIR + "/pre/{sample}/NanoStats.txt",
        pre_ybl=QC_DIR + "/pre/{sample}/Yield_By_Length.png",     
        pre_lvq=QC_DIR + "/pre/{sample}/LengthvsQualityScatterPlot_dot.png"  
    conda: 
        "../envs/01_qc.yml"
    threads: THREADS 
    shell: 
        """
        mkdir -p {QC_DIR}/pre/{wildcards.sample}
        NanoPlot --fastq {input.fastq} -o {QC_DIR}/pre/{wildcards.sample} --plots dot --format png --dpi 200 -t {threads} 
        """
# ------------ FILTRADO --------------------------------------------------------------------------
rule fitlong:
    input:
        fastq=RAW_DIR + "/combined/{sample}.fastq.gz"
    output:
        fastq_filt=FILT_DIR + "/{sample}_post_qc.fastq.gz",
    conda: 
        "../envs/01_qc.yml"	
    shell:
        """
        mkdir -p {FILT_DIR}
        filtlong --min_length {MIN_LEN} --min_mean_q {MIN_Q} {input.fastq} | gzip > {output.fastq_filt}
        """
# ------------- QC POST ----------------------------------------------------------------
rule nanoplot_post:
    input:
        fastq=FILT_DIR + "/{sample}_post_qc.fastq.gz"
    output:
        post_stats=QC_DIR + "/post/{sample}/NanoStats.txt",
        post_ybl=QC_DIR + "/post/{sample}/Yield_By_Length.png",     
        post_lvq=QC_DIR + "/post/{sample}/LengthvsQualityScatterPlot_dot.png" 
    conda: 
        "../envs/01_qc.yml"
    threads: THREADS
    shell:
        """
        mkdir -p {QC_DIR}/post/{wildcards.sample}
        NanoPlot --fastq {input.fastq} -o {QC_DIR}/post/{wildcards.sample} --plots dot --format png --dpi 200 -t {threads}
        """
# ---------------- REPORTE --------------------------------------------------------------------------
rule render_qc_report:
    input:
        pre_stats=expand(QC_DIR + "/pre/{sample}/NanoStats.txt", sample=SAMPLES),
        post_stats=expand(QC_DIR + "/post/{sample}/NanoStats.txt", sample=SAMPLES),
        pre_ybl=expand(QC_DIR + "/pre/{sample}/Yield_By_Length.png", sample=SAMPLES),
        pre_lvq=expand(QC_DIR + "/pre/{sample}/LengthvsQualityScatterPlot_dot.png", sample=SAMPLES),
        post_ybl=expand(QC_DIR + "/post/{sample}/Yield_By_Length.png", sample=SAMPLES),
        post_lvq=expand(QC_DIR + "/post/{sample}/LengthvsQualityScatterPlot_dot.png", sample=SAMPLES)
    output:
        html=HTML_DIR + "/01_qc.html"
    conda: 
        "../envs/render.yml"
    shell:
       """
        mkdir -p {HTML_DIR}
        Rscript --vanilla -e '
        rmarkdown::render(
          "analysis/01_qc.Rmd",
          params = list(
            pre_stats = "{input.pre_stats}",
            post_stats = "{input.post_stats}"
          ),
          knit_root_dir = getwd()
        )
        '
        mv analysis/01_qc.html {HTML_DIR}/01_qc.html
        """
# ---------- REPORTE PDF ------------------------------------------------------------------------------
rule html_to_pdf:
    input:
        html=HTML_DIR + "/01_qc.html"
    output:
        pdf=HTML_DIR + "/01_qc.pdf"
    conda:
        "../envs/render.yml"
    shell:
        """
        weasyprint {input.html} {output.pdf}
        """
