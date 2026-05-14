#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy import sparse


def load_reference(ref_h5ad: Path) -> ad.AnnData:
    return ad.read_h5ad(ref_h5ad, backed="r")


def sample_code(sample_name: str) -> str:
    return sample_name.split("_", 1)[0]


def build_barcode_map(ref_h5ad: Path, sample_name: str) -> tuple[dict[str, str], list[str], pd.DataFrame]:
    adata = load_reference(ref_h5ad)
    obs = adata.obs.copy()
    obs_names = adata.obs_names.to_list()
    mask = obs["sample"].astype(str) == sample_name
    cells = [name for name, keep in zip(obs_names, mask) if keep]
    if not cells:
        raise ValueError(f"No cells found for sample {sample_name}")

    code = sample_code(sample_name)
    barcode_to_cell = {}
    for cell in cells:
        barcode = cell.rsplit("-", 1)[0] + "-1"
        barcode_to_cell[barcode] = cell

    sample_obs = obs.loc[cells].copy()
    return barcode_to_cell, cells, sample_obs


def emit_triplets(args: argparse.Namespace) -> int:
    barcode_to_cell, _, _ = build_barcode_map(args.reference_h5ad, args.sample)
    whitelist = set(barcode_to_cell)

    bam = pysam.AlignmentFile(args.bam_url, "rb")
    processed = 0
    kept = 0

    try:
        for rec in bam.fetch(until_eof=True):
            processed += 1

            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            if not rec.has_tag("CB") or not rec.has_tag("UB"):
                continue
            if not rec.has_tag("GN") and not rec.has_tag("GX"):
                continue
            if rec.has_tag("RE") and rec.get_tag("RE") != "E":
                continue

            barcode = rec.get_tag("CB")
            if barcode not in whitelist:
                continue

            gene = rec.get_tag("GN") if rec.has_tag("GN") else rec.get_tag("GX")
            if not gene or ";" in gene:
                continue

            cell_name = barcode_to_cell[barcode]
            umi = rec.get_tag("UB")
            try:
                sys.stdout.write(f"{cell_name}\t{gene}\t{umi}\n")
            except BrokenPipeError:
                return 0
            kept += 1

            if processed % args.progress_every == 0:
                print(
                    f"[emit-triplets] processed={processed:,} kept={kept:,}",
                    file=sys.stderr,
                    flush=True,
                )
    finally:
        try:
            bam.close()
        except ConnectionResetError:
            # Remote S3-backed BAM streams sometimes reset on close after the useful data was already read.
            print("[emit-triplets] warning: remote BAM close reset; ignoring close-time failure", file=sys.stderr, flush=True)

    print(f"[emit-triplets] processed={processed:,} kept={kept:,}", file=sys.stderr, flush=True)
    return 0


def counts_to_h5ad(args: argparse.Namespace) -> int:
    _, cells, sample_obs = build_barcode_map(args.reference_h5ad, args.sample)
    ref = load_reference(args.reference_h5ad)
    genes = ref.var_names.to_list()
    var = ref.var.copy()

    cell_index = {cell: i for i, cell in enumerate(cells)}
    gene_index = {gene: i for i, gene in enumerate(genes)}

    rows = []
    cols = []
    data = []

    with args.counts_tsv.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            cell, gene, count = line.split("\t")
            if cell not in cell_index or gene not in gene_index:
                continue
            rows.append(cell_index[cell])
            cols.append(gene_index[gene])
            data.append(int(count))

    matrix = sparse.csr_matrix((data, (rows, cols)), shape=(len(cells), len(genes)), dtype=np.int32)
    out = ad.AnnData(X=matrix, obs=sample_obs, var=var)
    out.write_h5ad(args.output_h5ad)
    print(
        f"[counts-to-h5ad] wrote {args.output_h5ad} with {out.n_obs} cells x {out.n_vars} genes and {matrix.nnz} nonzeros"
    )
    return 0


def merge_h5ad(args: argparse.Namespace) -> int:
    adatas = [ad.read_h5ad(path) for path in args.input_h5ads]
    merged = ad.concat(adatas, axis=0, join="same", merge="same")
    merged.write_h5ad(args.output_h5ad)
    print(f"[merge-h5ad] wrote {args.output_h5ad} with {merged.n_obs} cells x {merged.n_vars} genes")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Recover raw UMI counts from public 10x BAM files.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    emit = subparsers.add_parser("emit-triplets", help="Emit cell/gene/UMI triplets for one sample.")
    emit.add_argument("--reference-h5ad", type=Path, required=True)
    emit.add_argument("--sample", required=True)
    emit.add_argument("--bam-url", required=True)
    emit.add_argument("--progress-every", type=int, default=5_000_000)
    emit.set_defaults(func=emit_triplets)

    counts = subparsers.add_parser("counts-to-h5ad", help="Build one sample raw-count h5ad from counted UMIs.")
    counts.add_argument("--reference-h5ad", type=Path, required=True)
    counts.add_argument("--sample", required=True)
    counts.add_argument("--counts-tsv", type=Path, required=True)
    counts.add_argument("--output-h5ad", type=Path, required=True)
    counts.set_defaults(func=counts_to_h5ad)

    merge = subparsers.add_parser("merge-h5ad", help="Merge per-sample raw-count h5ads.")
    merge.add_argument("--input-h5ads", type=Path, nargs="+", required=True)
    merge.add_argument("--output-h5ad", type=Path, required=True)
    merge.set_defaults(func=merge_h5ad)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
