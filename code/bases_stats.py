import pandas as pd
import os
import re

nanostats_post = snakemake.input.nanostats_post
nanostats_post_host = snakemake.input.nanostats_post_host
out_csv = snakemake.output.csv


def to_float(x: str) -> float:
    x = x.strip().replace(",", "").replace("%", "")
    return float(x)


def parse_nanostats(path: str) -> dict:
    """
    Extrae al menos:
      - Number of reads
      - Total bases
    de NanoStats.txt (NanoPlot).
    """
    out = {"reads": None, "bases": None}

    with open(path) as fp:
        for line in fp:
            if "Number of reads" in line:
                out["reads"] = int(to_float(line.split(":", 1)[1]))
            elif "Total bases" in line:
                out["bases"] = int(to_float(line.split(":", 1)[1]))

    if out["reads"] is None or out["bases"] is None:
        raise ValueError(f"No pude parsear reads/bases desde: {path}")

    return out


# Map sample -> filepath
def sample_from_nanostats_path(p: str) -> str:
    # .../post/{sample}/NanoStats.txt  o  .../post_host/{sample}/NanoStats.txt
    return os.path.basename(os.path.dirname(p))


post_map = {sample_from_nanostats_path(p): p for p in nanostats_post}
post_host_map = {sample_from_nanostats_path(p): p for p in nanostats_post_host}

samples = sorted(set(post_map.keys()) & set(post_host_map.keys()))
if not samples:
    raise ValueError("No hay samples en común entre nanostats_post y nanostats_post_host.")

rows = []
for s in samples:
    post = parse_nanostats(post_map[s])
    dehost = parse_nanostats(post_host_map[s])

    removed_bases = post["bases"] - dehost["bases"]
    removed_reads = post["reads"] - dehost["reads"]

    pct_bases_removed = (removed_bases / post["bases"] * 100) if post["bases"] > 0 else 0.0
    pct_reads_removed = (removed_reads / post["reads"] * 100) if post["reads"] > 0 else 0.0

    rows.append({
        "sample": s,
        "post_reads": post["reads"],
        "post_bases": post["bases"],
        "dehost_reads": dehost["reads"],
        "dehost_bases": dehost["bases"],
        "removed_reads": removed_reads,
        "removed_bases": removed_bases,
        "pct_reads_removed": pct_reads_removed,
        "pct_bases_removed": pct_bases_removed,
    })

df = pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)
df.to_csv(out_csv, index=False)
