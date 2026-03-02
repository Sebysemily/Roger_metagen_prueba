import pandas as pd
import re
import os

flagstats = snakemake.input.flagstats
nanostats = snakemake.input.nanostats
output_csv = snakemake.output.csv


def to_float(x: str) -> float:
    x = x.strip().replace(",", "").replace("%", "")
    return float(x)


def to_int(x: str) -> int:
    return int(to_float(x))


def parse_sample_host(basename: str):
    """
    Acepta:
      - A1.host        -> sample=A1, host=host
      - A1_mutuki      -> sample=A1, host=mutuki
      - A1_palythoa_v1 -> sample=A1_palythoa, host=v1 (si solo quieres split final por _)
    """
    # preferir formato con punto: sample.host
    if "." in basename:
        parts = basename.split(".", 1)
        if len(parts) == 2 and parts[0] and parts[1]:
            return parts[0], parts[1]

    # fallback: formato con underscore: sample_host (split por el último _)
    if "_" in basename:
        parts = basename.rsplit("_", 1)
        if len(parts) == 2 and parts[0] and parts[1]:
            return parts[0], parts[1]

    raise ValueError(f"No puedo parsear sample/host desde basename='{basename}'")


# --------- Leer total reads desde NanoStats ----------
total_reads = {}
missing_totals = []

for f in nanostats:
    sample = os.path.basename(os.path.dirname(f))  # .../post/{sample}/NanoStats.txt
    found = False
    with open(f) as fp:
        for line in fp:
            if "Number of reads" in line:
                total_reads[sample] = to_int(line.split(":", 1)[1])
                found = True
                break
    if not found:
        missing_totals.append(f)

if missing_totals:
    raise ValueError(
        "No se encontró la línea 'Number of reads' en:\n" + "\n".join(missing_totals)
    )


# --------- Parsear flagstat correctamente ----------
data = []
bad_flagstats = []

# regex para la línea "mapped"
mapped_re = re.compile(r"([\d,]+)\s+\+\s+[\d,]+\s+mapped.*?\(([\d.]+)%")

for f in flagstats:
    basename = os.path.basename(f).replace(".flagstat", "")
    sample, host = parse_sample_host(basename)

    mapped = 0
    pct_flagstat = 0.0
    found_mapped_line = False

    with open(f) as fp:
        for line in fp:
            if " mapped" in line:
                m = mapped_re.search(line)
                if m:
                    mapped = to_int(m.group(1))
                    pct_flagstat = to_float(m.group(2))
                    found_mapped_line = True
                    break

    if not found_mapped_line:
        bad_flagstats.append(f)

    total = total_reads.get(sample, 0)
    pct_real = (mapped / total * 100) if total > 0 else 0.0

    data.append([sample, host, mapped, total, pct_flagstat, pct_real])

if bad_flagstats:
    raise ValueError(
        "No pude encontrar/parsear la línea 'mapped' en estos flagstat:\n"
        + "\n".join(bad_flagstats)
    )

df = pd.DataFrame(
    data,
    columns=["sample", "host", "mapped_reads", "total_reads", "pct_flagstat", "pct_real"]
).sort_values(["sample", "host"]).reset_index(drop=True)

df.to_csv(output_csv, index=False)
