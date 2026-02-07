# Inputs:
# - Raw sequencing data in Data/Raw_seq/..fastq.gz
# Outputs:
# - Filtered fastq in Data/fitlon/...fastq.gz
 
 rule nanoplot_pre:
    input:
	"data/raw_seq/{sample}.fast.gz"
    output:
	pre_stats="results/outputs/qc/pre/{sample}/NanoStats.txt",
	pre_ybl="results/outputs/qc/pre/{sample}/Yield_By_Length.html",
	pre_lvq="results/outputs/qc/pre/{sample}/LengthvsQualityScatterPlot_kde.html"
    conda: 
	"envs/01_qc/nanoplot.yml" 
    threads: config["threads"]
    shell: 
	"""
	mkdir -p results/qc/pre/{wildcards.sample}
	NanoPlot --fastq {input} -o results/outputs/qc/pre/{wildcards.sample} --threads {threads} 
	"""
	
