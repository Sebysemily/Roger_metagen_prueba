# Inputs:
# - Raw sequencing data in Data/Raw_seq/..fastq.gzip
# Outputs:
# - Filtered fastq in Data/fitlon/...fastq.gzip
 
 rule nanopot_pre:
    input:
	fastq = "data/{id}.fast.gzip"
    output:
	stats="results/outputs/nanoplot_pre_qc/NanoStats.txt",
	Yield_By_Length_png="results/outputs/nanoplot_pre_qc/Yield_By_Length.png",
	Yield_By_Length_html="results/outputs/nanoplot_pre_qc/Yield_By_Length.html",
	W_H_Length_png="results/outputs/nanoplot_pre_qc/WeightedHistogramReadlength.png",
	W_H_Length_html="results/outputs/nanoplot_pre_qc/WeightedHistogramReadle.html",
	L_vs_Q_png="results/outptus/nanoplot_pre_qc/LengthvsQualityScatterPlot_kde.png",
	L_vs_Q_html="results/outptus/nanoplot_pre_qc/LengthvsQualityScatterPlot.html"
    conda: 
	"envs/01_qc/nanoplot.yml" 
    shell: 
	"""
	Nanoplot --fastq {input} -o results/outputs/nanoplot_pre_qc/{wildcards.sample} 
	"""
	
