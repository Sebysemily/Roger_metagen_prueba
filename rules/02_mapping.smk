import pandas as pd

configfile: "config/config.yml"

SAMPLES = pd.read_csv("config/samples.tsv", sep="\t")["sample"].tolist()

FILT_DIR = config["paths"]["filtered_dir"]
HTML_DIR = config["paths"]["html_dir"]
REF_DIR = config["paths"]["ref_dir"]
MAPPING_DIR = config["paths"]["mapping_dir"]
QC_DIR= config["paths"]["qc_dir"]
DEHOSTED_DIR= config["pahts"]["dehosted_dir"]

THREADS = config["mapping"]["threads"]
EXTRA = config["mapping"]["extra_options"]
HOSTS = config["mapping"]["hosts"]

#-----------Indexar genomas host------------------------------------------------
rule index_host_genome:
    input: 
        genome=REF_DIR +"/{host}.fna"
    output:
        index= REF_DIR +"/{host}.mmi"
    conda:
        "../envs/mapping.yml"
    threads: THREADS
    shell:
        """
        minimap2 -x mpa-ont -d {output.index} {input.genome}
        """
#----------- Mappear a host ----------------------------------------------------
rule map_sample_to_host:
    input:
        reads=FILT_DIR + "/{sample}_post_qc.fastq.gz",
        index=REF_DIR +"/{host}.mmi"
    output:
        bam=MAPPING_DIR + "/bam/{sample}_{host}.bam"
    conda:
        "../envs/02_mapping.yml"
    threads: THREADS
    shell:
        """
        minimap2 -t {threads} -ax map-ont {EXTRA} {input.index} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """
#----------Calcular estadisticas de mappeo--------------------------------------
rule flagstat:
    input:
        bam=MAPPING_DIR + "/bam/{sample}_{host}.bam"
    output:
        stats=MAPPING_DIR + "/flagstat/{sample}_{host}.flagstat"
    conda:
        "../envs/02_mapping.yml"
    shell:
        "samtools flagstat {input.bam} > {output.stats}"
#-------Generar CSV con % de mapeo----------------------------------------------
rule mapping_stats:
    input: 
        flagstats=expand(MAPPING_DIR + "/flagstat/{sample}_{host}.flagstat",
                         sample=SAMPLES, host=HOSTS),
        nanostats=expand(QC_DIR + "/post/{sample}/NanoStats.txt", sample=SAMPLES)
    output:
        csv=MAPPING_DIR + "/stats_csv/mapping_stats.csv"
    script:
        "code/gather_mapping_stats.py"
#------Extraer secuencias mappeadas---------------------------------------------
rule extract_mapped_ids:
    input:
        bam=MAPPING_DIR + "/bam/{sample}_{host}.bam"
    output:
        ids=MAPPING_DIR + "/mapped_ids/{sample}_{host}.txt"
    conda:
        "../envs/02_mapping.yml"
    shell:
        """
        samtools view -F 4 {input.bam} | cut -f1 | sort -u > {output.ids}
        """
#------------Unir seq mapeadas--------------------------------------------------
rule union_mapped_ids:
    input:
        ids=expand(MAPPING_DIR + "/mapped_ids/{{sample}}_{host}.txt",
                   host=HOSTS)
    output:
        blacklist=MAPPING_DIR + "/mapped_ids/{sample}_all_hosts.txt"
    run:
        all_ids = set()
        for f in input.ids:
            with open(f) as fp:
                all_ids.update(line.strip() for line in fp)
        with open(output.blacklist, "w") as out:
            for rid in sorted(all_ids):
                out.write(rid + "\n")
#--------- quitar host----------------------------------------------------------
rule dehost_fastq:
    input:
        fastq=FILT_DIR + "/{sample}_post_qc.fastq.gz",
        blacklist=MAPPING_DIR + "/mapped_ids/{sample}_all_hosts.txt"
    output:
        clean=DEHOSTED_DIR + "/{sample}_dehosted.fastq.gz"
    conda:
        "../envs/seqkit.yml"   
    shell:
        """
        seqkit grep -v -f {input.blacklist} {input.fastq} | gzip > {output.clean}
        """

