include: "rules/01_qc_local.smk"

rule all:
    input:
      "results/results_docs/01_qc.html"
