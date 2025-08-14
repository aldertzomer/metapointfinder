#!/usr/bin/env python3
import sys
import re
import random
import argparse
from typing import Dict, List, Tuple, Optional

# ---------------- FASTA ----------------
def read_fasta_by_accession(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    acc = None
    buf: List[str] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf).replace(" ", "").upper()
                buf = []
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 2 else parts[0].split()[0]
            else:
                buf.append(line.strip())
        if acc is not None:
            seqs[acc] = "".join(buf).replace(" ", "").upper()
    return seqs

# --------------- Tokens ----------------
MUT_RE = re.compile(r"^([A-Z]+)(\d+)([A-Z\-]+)$")

def parse_token(token: str) -> Optional[Tuple[str, int, str]]:
    t = token.strip().upper()
    m = MUT_RE.match(t)
    if not m:
        return None
    ref_seg, pos_s, alt_seg = m.group(1), m.group(2), m.group(3)
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
    s = seq
    L = len(s)
    ref_len = len(ref_seg)
    idx0 = pos - 1

    def matches_at(i: int) -> bool:
        return 0 <= i <= L - ref_len and s[i:i+ref_len] == ref_seg

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
        return s[:idx0] + s[idx0 + ref_len:]
    if ref_len == 1 and alt_seg.startswith(ref_seg) and len(alt_seg) > 1:
        ins = alt_seg[1:]
        return s[:idx0 + 1] + ins + s[idx0 + 1:]
    return s[:idx0] + alt_seg + s[idx0 + ref_len:]

# --------------- IO helpers ------------
def wrap60(s: str):
    for i in range(0, len(s), 60):
        yield s[i:i+60]

# --------------- Main ------------------
def main():
    ap = argparse.ArgumentParser(description="Generate wildtype + mutant protein sequences from AMR mutation tokens.")
    ap.add_argument("tsv", help="AMRProt-mutation_underscore_combined.tsv")
    ap.add_argument("fasta", help="AMRProt_mutation.fa")
    ap.add_argument("out", help="Output FASTA file")
    ap.add_argument("--per-row", type=int, default=1,
                    help="Number of mutant+WT pairs to generate per TSV line (default: 1)")
    ap.add_argument("--with-replacement", action="store_true",
                    help="Sample mutation tokens with replacement (default: without)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    acc2seq = read_fasta_by_accession(args.fasta)
    if not acc2seq:
        sys.stderr.write(f"[ERROR] No sequences parsed from {args.fasta}\n")
        sys.exit(1)

    header_id = 1
    out_records: List[Tuple[str, str]] = []

    with open(args.tsv, "r") as fh:
        header = fh.readline()
        if not header:
            sys.stderr.write(f"[ERROR] Empty TSV: {args.tsv}\n")
            sys.exit(1)

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            cls, acc, changes_str = parts[0], parts[1], parts[2]
            if acc not in acc2seq:
                sys.stderr.write(f"[WARN] Accession {acc} not found in FASTA, skipping.\n")
                continue

            wt_seq = acc2seq[acc]
            mut_tokens = [t for t in changes_str.split(",") if parse_token(t)]

            if not mut_tokens:
                sys.stderr.write(f"[WARN] No valid mutations for {acc}, skipping.\n")
                continue

            sample_size = min(args.per_row, len(mut_tokens)) if not args.with_replacement else args.per_row
            chosen_tokens = rng.sample(mut_tokens, sample_size) if not args.with_replacement else \
                            [rng.choice(mut_tokens) for _ in range(sample_size)]

            for tok in chosen_tokens:
                ref_seg, pos, alt_seg = parse_token(tok)
                try:
                    mut_seq = apply_mutation(wt_seq, ref_seg, pos, alt_seg)
                except ValueError as e:
                    sys.stderr.write(f"[WARN] Could not apply {tok} to {acc}: {e}\n")
                    continue

                out_records.append((f">{acc}|{cls}|RES={tok}|{header_id}", mut_seq))
                header_id += 1
                out_records.append((f">{acc}|{cls}|WT|{header_id}", wt_seq))
                header_id += 1

    with open(args.out, "w") as out_fh:
        for hdr, seq in out_records:
            out_fh.write(hdr + "\n")
            out_fh.write("\n".join(wrap60(seq)) + "\n")

if __name__ == "__main__":
    main()
