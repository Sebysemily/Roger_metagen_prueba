import pandas as pd
import re
import os

# Los objetos input, output, etc. están disponibles en snakemake
flagstats = snakemake.input.flagstats   # lista
nanostats = snakemake.input.nanostats   # lista
output_csv = snakemake.output.csv

# Mismo código que antes, pero usando estas variables directamente
total_reads = {}
for f in nanostats:
    sample = os.path.basename(os.path.dirname(f))
    with open(f) as fp:
        for line in fp:
            if "Number of reads" in line:
                total = int(line.split(":")[1].strip())
                total_reads[sample] = total
                break

data = []
for f in flagstats:
    basename = os.path.basename(f).replace(".flagstat", "")
    parts = basename.split("_")
    sample = "_".join(parts[:-1])
    host = parts[-1]
    with open(f) as fp:
        first_line = fp.readline()
        match = re.search(r"(\d+)\s+\+\s+\d+\s+mapped.*?\((\d+\.?\d*)%\)", first_line)
        if match:
            mapped = int(match.group(1))
            pct_flagstat = float(match.group(2))
        else:
            mapped = 0
            pct_flagstat = 0.0
    total = total_reads.get(sample, 0)
    pct_real = (mapped / total * 100) if total > 0 else 0
    data.append([sample, host, mapped, total, pct_flagstat, pct_real])

df = pd.DataFrame(data, columns=["sample", "host", "mapped_reads",
                                 "total_reads", "pct_flagstat", "pct_real"])
df.to_csv(output_csv, index=False)
