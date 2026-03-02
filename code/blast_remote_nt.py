#!/usr/bin/env python3
"""
Snakemake script: blast_remote_nt.py

- Divide un FASTA en chunks
- Ejecuta blastn -remote contra NCBI (nt por defecto)
- Maneja rate-limit / fallos de conexión con reintentos + backoff
- Puede "resume": si ya existe el .nt.tsv de un chunk y tiene contenido, lo salta
- Puede "fail_soft": si un chunk falla tras reintentos, crea un .nt.tsv vacío y continúa

Outputs:
  - {output.tsv}: concatenación de todos los chunks
  - {output.summary}: resumen simple por taxid (o "NA" si no hay taxid)
"""

import csv
import random
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import List, Tuple

# Snakemake injects `snakemake`
# noqa: F821


def log(msg: str) -> None:
    print(f"[blast_remote] {msg}", flush=True)


def run_cmd(cmd: List[str]) -> Tuple[int, str, str]:
    """Run command, capture stdout/stderr."""
    p = subprocess.run(cmd, text=True, capture_output=True)
    return p.returncode, p.stdout, p.stderr


def is_retryable_blast_error(stderr: str) -> bool:
    s = stderr.lower()
    patterns = [
        "bad_request",
        "could not queue request",
        "request submission failed",
        "connection stream is in bad state",
        "blast4-request",
        "service unavailable",
        "too many requests",
        "timeout",
        "temporarily",
        "try again",
    ]
    return any(p in s for p in patterns)


def split_fasta(in_fa: Path, outdir: Path, chunk_size: int) -> List[Path]:
    """Split FASTA into chunk files with chunk_size sequences each."""
    outdir.mkdir(parents=True, exist_ok=True)
    chunks: List[Path] = []

    def new_chunk(idx: int) -> Path:
        return outdir / f"{in_fa.stem}.part_{idx:03d}.fasta"

    idx = 1
    seqs_in_chunk = 0
    out_handle = None

    def open_new():
        nonlocal idx, seqs_in_chunk, out_handle
        if out_handle:
            out_handle.close()
        chunk_path = new_chunk(idx)
        chunks.append(chunk_path)
        out_handle = chunk_path.open("w")
        seqs_in_chunk = 0
        idx += 1

    open_new()

    with in_fa.open() as f:
        header = None
        seq_lines = []
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    out_handle.write(header)
                    out_handle.writelines(seq_lines)
                    seqs_in_chunk += 1
                    if seqs_in_chunk >= chunk_size:
                        open_new()
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            out_handle.write(header)
            out_handle.writelines(seq_lines)

    if out_handle:
        out_handle.close()

    chunks = [c for c in chunks if c.stat().st_size > 0]
    return chunks


def blast_one_chunk(
    chunk_fa: Path,
    out_tsv: Path,
    db: str,
    task: str,
    evalue: str,
    max_target_seqs: int,
    include_scinames: bool,
    retries: int,
    sleep_s: float,
    jitter_s: float,
    backoff: float,
    resume: bool,
    fail_soft: bool,
) -> None:
    if resume and out_tsv.exists() and out_tsv.stat().st_size > 0:
        log(f"resume: skip existing {out_tsv.name}")
        return

    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    if include_scinames:
        outfmt = "6 qseqid sacc staxids sscinames pident length qcovs evalue bitscore stitle"
    else:
        outfmt = "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle"

    cmd = [
        "blastn",
        "-remote",
        "-db",
        db,
        "-task",
        task,
        "-query",
        str(chunk_fa),
        "-max_target_seqs",
        str(max_target_seqs),
        "-max_hsps",
        "1",
        "-evalue",
        str(evalue),
        "-outfmt",
        outfmt,
        "-out",
        str(out_tsv),
    ]

    attempt = 0
    delay = sleep_s

    while True:
        attempt += 1
        rc, _stdout, stderr = run_cmd(cmd)

        if rc == 0:
            return

        retryable = is_retryable_blast_error(stderr)
        log(f"chunk {chunk_fa.name} failed (rc={rc}) attempt {attempt}/{retries}. retryable={retryable}")

        tail = "\n".join(stderr.strip().splitlines()[-3:])
        if tail:
            log(f"stderr tail:\n{tail}")

        if (not retryable) or (attempt >= retries):
            if fail_soft:
                log(f"fail_soft: writing empty result for {chunk_fa.name}")
                out_tsv.write_text("")
                return
            raise RuntimeError(f"blastn failed for {chunk_fa} after {attempt} attempts (rc={rc}). stderr:\n{stderr}")

        sleep_for = delay + random.uniform(0, jitter_s)
        log(f"sleep {sleep_for:.1f}s then retry...")
        time.sleep(sleep_for)
        delay *= backoff


def concat_tsv(chunk_tsvs: List[Path], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as out:
        for p in chunk_tsvs:
            if p.exists():
                out.write(p.read_text())


def summarize(out_tsv: Path, out_csv: Path) -> None:
    total_hits = 0
    q_with_hit = set()
    tax_counter = Counter()

    with out_tsv.open() as f:
        for line in f:
            if not line.strip():
                continue
            total_hits += 1
            cols = line.rstrip("\n").split("\t")
            qseqid = cols[0]
            staxids = cols[2] if len(cols) > 2 else "NA"
            q_with_hit.add(qseqid)

            if staxids.strip() == "" or staxids == "N/A":
                tax_counter["NA"] += 1
            else:
                tax_counter[staxids.split(";")[0]] += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["total_hits", total_hits])
        w.writerow(["unique_queries_with_hit", len(q_with_hit)])
        w.writerow([])
        w.writerow(["top_taxid", "hit_count"])
        for taxid, cnt in tax_counter.most_common(20):
            w.writerow([taxid, cnt])


def main() -> None:
    in_fasta = Path(snakemake.input["fasta"])       # noqa: F821
    out_tsv = Path(snakemake.output["tsv"])         # noqa: F821
    out_summary = Path(snakemake.output["summary"]) # noqa: F821

    p = snakemake.params  # noqa: F821
    chunks_dir = Path(p["chunks_dir"])

    db = str(p["db"])
    task = str(p["task"])
    evalue = str(p["evalue"])
    max_target_seqs = int(p["max_target_seqs"])
    chunk_size = int(p["chunk_size"])

    include_scinames = bool(p.get("include_scinames", False))
    sleep_s = float(p.get("sleep_s", 8))
    jitter_s = float(p.get("jitter_s", 5))
    retries = int(p.get("retries", 8))
    backoff = float(p.get("backoff", 1.8))
    fail_soft = bool(p.get("fail_soft", True))
    resume = bool(p.get("resume", True))

    chunks_dir.mkdir(parents=True, exist_ok=True)

    existing = sorted(chunks_dir.glob(f"{in_fasta.stem}.part_*.fasta"))
    if existing:
        chunks = existing
        log(f"using existing {len(chunks)} chunk FASTAs from {chunks_dir}")
    else:
        log(f"split {in_fasta} into chunks of {chunk_size}")
        chunks = split_fasta(in_fasta, chunks_dir, chunk_size)
        log(f"created {len(chunks)} chunk FASTAs")

    chunk_tsvs: List[Path] = []
    for ch in chunks:
        ch_tsv = ch.with_suffix(ch.suffix + ".nt.tsv")  # *.fasta.nt.tsv
        chunk_tsvs.append(ch_tsv)

        time.sleep(sleep_s + random.uniform(0, jitter_s))

        blast_one_chunk(
            chunk_fa=ch,
            out_tsv=ch_tsv,
            db=db,
            task=task,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            include_scinames=include_scinames,
            retries=retries,
            sleep_s=sleep_s,
            jitter_s=jitter_s,
            backoff=backoff,
            resume=resume,
            fail_soft=fail_soft,
        )

    concat_tsv(chunk_tsvs, out_tsv)
    summarize(out_tsv, out_summary)
    log(f"done: {out_tsv} ({out_tsv.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
