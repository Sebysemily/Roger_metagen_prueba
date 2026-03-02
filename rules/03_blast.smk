import pandas as pd

configfile: "config/config.yml"

SAMPLES = pd.read_csv("config/samples.tsv", sep="\t")["sample"].tolist()

# -------------------------
# Paths (robusto a KeyError)
# -------------------------

REF_DIR      = config["paths"]["ref_dir"]
DEHOSTED_DIR = config["paths"]["dehosted_dir"]
BLAST_DIR    = config["paths"]["blast_dir"]

THREADS  = config["mapping"]["threads"]
HOSTS    = config["mapping"]["hosts"]

HOST_REF_FASTA = REF_DIR + "/hosts_combined.fna"
HOST_BLASTDB   = REF_DIR + "/hosts_combined_blastdb"

# -------------------------
# Local BLAST (removed vs hosts)
LOCAL_TASK       = config["blast"].get("local_task", "megablast")
LOCAL_EVALUE     = str(config["blast"].get("local_evalue", "1e-20"))
LOCAL_MAX_TARGET = int(config["blast"].get("local_max_target_seqs", 25))

# Remote BLAST params
REMOTE_DB          = config["blast"].get("remote_db", "nt")
REMOTE_TASK        = config["blast"].get("remote_task", "megablast")
REMOTE_EVALUE      = str(config["blast"].get("remote_evalue", "1e-10"))
REMOTE_MAX_TARGET  = int(config["blast"].get("remote_max_target_seqs", 10))
REMOTE_MAX_READS   = int(config["blast"].get("remote_max_reads", 200))
REMOTE_CHUNK_SIZE  = int(config["blast"].get("remote_chunk_size", 20))



# -------------------------------------------------------------------------
# 1) Crear BLAST DB local de hosts
# -------------------------------------------------------------------------
rule make_hosts_blastdb:
    input:
        ref=HOST_REF_FASTA
    output:
        # makeblastdb crea varios archivos; usamos .nin como sentinel
        nin=HOST_BLASTDB + ".nin"
    conda:
        "../envs/03_blast.yml"
    shell:
        r"""
        makeblastdb -in {input.ref} -dbtype nucl -out {HOST_BLASTDB}
        """

# -------------------------------------------------------------------------
# 2) Convertir removed.fastq.gz a fasta (para blast local)
# -------------------------------------------------------------------------
rule removed_to_fasta:
    input:
        fq=os.path.join(DEHOSTED_DIR, "{sample}_removed.fastq.gz")
    output:
        fa=os.path.join(BLAST_DIR, "tmp", "{sample}.removed.fasta")
    conda:
        "../envs/03_blast.yml"
    shell:
        r"""
        mkdir -p {BLAST_DIR}/tmp
        seqkit fq2fa {input.fq} > {output.fa}
        """

# -------------------------------------------------------------------------
# 3) BLAST local: removed vs hosts_combined_blastdb
# -------------------------------------------------------------------------
rule blast_removed_vs_hosts_local:
    input:
        fa=os.path.join(BLAST_DIR, "tmp", "{sample}.removed.fasta"),
        nin=HOST_BLASTDB + ".nin"
    output:
        tsv=os.path.join(BLAST_DIR, "local_hosts", "{sample}.removed.vs_hosts.tsv")
    conda:
        "../envs/03_blast.yml"
    threads: config["mapping"].get("threads", 10)
    shell:
        r"""
        mkdir -p {BLAST_DIR}/local_hosts
        blastn \
          -query {input.fa} \
          -db {HOST_BLASTDB} \
          -task {LOCAL_TASK} \
          -evalue {LOCAL_EVALUE} \
          -max_target_seqs {LOCAL_MAX_TARGET} \
          -num_threads {threads} \
          -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
          > {output.tsv}
        """

# -------------------------------------------------------------------------
# 4) Preparar dehosted fasta subsample (para BLAST remoto nt)
# -------------------------------------------------------------------------
rule dehosted_to_fasta_subsample:
    input:
        fq=os.path.join(DEHOSTED_DIR, "{sample}_dehosted.fastq.gz")
    output:
        fa=os.path.join(BLAST_DIR, "tmp", "{sample}.dehosted.subset.fasta")
    conda:
        "../envs/03_blast.yml"
    shell:
        r"""
        mkdir -p {BLAST_DIR}/tmp
        seqkit sample -n {REMOTE_MAX_READS} {input.fq} | seqkit fq2fa > {output.fa}
        """

# -------------------------------------------------------------------------
# 5) BLAST remoto nt por chunks (robusto: retry/backoff + resume)
#    IMPORTANTE:
#      Para no saturar NCBI, corre Snakemake con:
#        snakemake ... --resources ncbi_remote=1
# -------------------------------------------------------------------------
rule blast_dehosted_vs_remote:
    input:
        fasta=os.path.join(BLAST_DIR, "tmp", "{sample}.dehosted.subset.fasta")
    output:
        tsv=os.path.join(BLAST_DIR, "remote_nt", "{sample}.dehosted.vs_remote.tsv"),
        summary=os.path.join(BLAST_DIR, "remote_nt", "{sample}.dehosted.vs_remote.summary.csv"),
    params:
        db=REMOTE_DB,
        task=REMOTE_TASK,
        evalue=REMOTE_EVALUE,
        max_target_seqs=REMOTE_MAX_TARGET,
        chunk_size=REMOTE_CHUNK_SIZE,
        chunks_dir=os.path.join(BLAST_DIR, "remote_nt", "chunks", "{sample}"),
        include_scinames=config["blast"].get("remote_include_scinames", False),
        sleep_s=float(config["blast"].get("remote_sleep_s", 8)),
        jitter_s=float(config["blast"].get("remote_jitter_s", 5)),
        retries=int(config["blast"].get("remote_retries", 8)),
        backoff=float(config["blast"].get("remote_backoff_factor", 1.8)),
        fail_soft=bool(config["blast"].get("remote_fail_soft", True)),
        resume=bool(config["blast"].get("remote_resume", True)),
    conda:
        "../envs/03_blast.yml"
    threads: 1
    resources:
        ncbi_remote=1
    script:
        "../code/blast_remote_nt.py"
