#!/usr/bin/env python3
"""
make_dna_mutants.py — Generate wildtype + mutant DNA sequences from AMR DNA tokens.

Input TSV (tab-delimited): AMR_DNA_underscore_combined.tsv
  Columns (4):
    class    accession    sequence    changes_str
  - `changes_str` contains comma-separated tokens like: A198-, G523T, AA101AT, etc.

Output:
  FASTA with per-row mutant + WT pairs:
    >ACCESSION|CLASS|RES=<token>|<id>
    <mutant sequence>
    >ACCESSION|CLASS|WT|<id>
    <wildtype sequence>

Usage:
  python make_dna_mutants.py AMR_DNA_underscore_combined.tsv out.fna \
      --per-row 3 --seed 42
  # sampling without replacement across tokens by default; use --with-replacement to allow repeats.
"""

import sys
import re
import random
import argparse
from typing import List, Tuple, Optional

# ---------------- Tokens ----------------
# Allow A/C/G/T/N blocks to be safe (some DBs include N). Hyphen is deletion.
MUT_RE = re.compile(r"^([ACGTN]+)(\d+)([ACGTN\-]+)$", re.IGNORECASE)

def parse_token(token: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse DNA mutation token e.g. A198-, G523T, AA101AT, A100AT (keep A, insert T).
    Returns (ref_seg, pos (1-based), alt_seg) uppercase, or None if invalid.
    """
    t = token.strip().upper()
    if not t:
        return None
    m = MUT_RE.match(t)
    if not m:
        return None
    ref_seg, pos_s, alt_seg = m.group(1).upper(), m.group(2), m.group(3).upper()
    try:
        pos = int(pos_s)
    except ValueError:
        return None
    return ref_seg, pos, alt_seg

def apply_mutation(seq: str,
                   ref_seg: str,
                   pos: int,
                   alt_seg: str,
                   search_window: int = 3) -> str:
    """
    Apply a DNA mutation token to 'seq'.
      - pos is 1-based index of the FIRST base of ref_seg (reference numbering).
      - If exact match not at pos-1, search within ±search_window for ref_seg.
    Supports:
      * deletion: alt_seg == "-"
      * substitution / equal-length replacement
      * simple insert-after: len(ref)==1 and alt startswith(ref) and len(alt)>1
    Returns mutated sequence. Raises ValueError if ref_seg not found nearby.
    """
    s = seq
    L = len(s)
    ref_len = len(ref_seg)
    idx0 = pos - 1

    if ref_len <= 0 or idx0 < 0 or idx0 >= L:
        raise ValueError(f"Bad indices for token at pos {pos} with ref_len={ref_len} on len={L}")

    def matches_at(i: int) -> bool:
        return 0 <= i <= L - ref_len and s[i:i+ref_len] == ref_seg

    # try exact, then a small ± window
    if not matches_at(idx0):
        found = False
        for delta in range(-search_window, search_window + 1):
            if delta == 0:
                continue
            j = idx0 + delta
            if matches_at(j):
                idx0 = j
                found = True
                break
        if not found:
            raise ValueError(f"Ref segment '{ref_seg}' not found near pos {pos} (len={L}).")

    if alt_seg == "-":
        # deletion of the ref_seg block
        return s[:idx0] + s[idx0 + ref_len:]

    if ref_len == 1 and alt_seg.startswith(ref_seg) and len(alt_seg) > 1:
        # insertion-after shorthand: e.g., A -> AT (insert T after A)
        ins = alt_seg[1:]
        return s[:idx0 + 1] + ins + s[idx0 + 1:]

    # general replacement (covers substitutions and multi-base replacements)
    return s[:idx0] + alt_seg + s[idx0 + ref_len:]

# --------------- IO helpers ------------
def wrap60(s: str):
    for i in range(0, len(s), 60):
        yield s[i:i+60]

def parse_dna_tsv(tsv_path: str):
    """
    Yield tuples of (cls, acc, seq, changes_str) from TSV.
    Expects 4 tab-separated columns with header on first line.
    """
    with open(tsv_path, "r") as fh:
        header = fh.readline()
        if not header:
            raise RuntimeError(f"Empty TSV: {tsv_path}")
        # Simple validation (optional)
        # columns = header.rstrip("\n").split("\t")
        # if len(columns) < 4: raise RuntimeError("TSV must have at least 4 columns")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            cls, acc, seq, changes_str = parts[0], parts[1], parts[2], parts[3]
            yield cls, acc, seq, changes_str

# --------------- Main ------------------
def main():
    ap = argparse.ArgumentParser(description="Generate wildtype + mutant DNA sequences from AMR DNA tokens.")
    ap.add_argument("tsv", help="AMR_DNA_underscore_combined.tsv")
    ap.add_argument("out", help="Output FASTA file")
    ap.add_argument("--per-row", type=int, default=1,
                    help="Number of mutant+WT pairs to generate per TSV line (default: 1)")
    ap.add_argument("--with-replacement", action="store_true",
                    help="Sample mutation tokens with replacement (default: without)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    ap.add_argument("--search-window", type=int, default=3,
                    help="± columns to search if ref segment not found at exact pos (default: 3)")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    header_id = 1
    out_records = []  # list of (header, seq)

    for cls, acc, seq, changes_str in parse_dna_tsv(args.tsv):
        wt_seq = seq.strip().upper().replace(" ", "")
        if not wt_seq:
            sys.stderr.write(f"[WARN] Empty sequence for accession {acc}, skipping.\n")
            continue

        tokens = [t.strip() for t in changes_str.split(",") if t.strip()]
        parsed_tokens = [t for t in tokens if parse_token(t) is not None]
        if not parsed_tokens:
            sys.stderr.write(f"[WARN] No valid mutations for {acc}, skipping.\n")
            continue

        # choose how many tokens to generate per row
        sample_size = args.per_row if args.with_replacement else min(args.per_row, len(parsed_tokens))
        chosen = [rng.choice(parsed_tokens) for _ in range(sample_size)] if args.with_replacement \
                 else rng.sample(parsed_tokens, sample_size)

        for tok in chosen:
            p = parse_token(tok)
            if p is None:
                continue
            ref_seg, pos, alt_seg = p
            try:
                mut_seq = apply_mutation(wt_seq, ref_seg, pos, alt_seg, search_window=args.search_window)
            except ValueError as e:
                sys.stderr.write(f"[WARN] Could not apply {tok} to {acc}: {e}\n")
                continue

            # mutant
            out_records.append((f">{acc}|{cls}|RES={tok}|{header_id}", mut_seq))
            header_id += 1
            # wildtype
            out_records.append((f">{acc}|{cls}|WT|{header_id}", wt_seq))
            header_id += 1

    if not out_records:
        sys.stderr.write("[INFO] No records written (no valid mutations applied?).\n")

    with open(args.out, "w") as out_fh:
        for hdr, seq in out_records:
            out_fh.write(hdr + "\n")
            out_fh.write("\n".join(wrap60(seq)) + "\n")

if __name__ == "__main__":
    main()
