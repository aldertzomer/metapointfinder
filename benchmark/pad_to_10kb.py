#!/usr/bin/env python3
"""
pad_to_10kb.py — Pad/truncate FASTA sequences to exactly 10,000 bp.

- If a sequence is shorter than 10,000, add random A/C/G/T bases on both sides.
- If a sequence is longer than 10,000, trim centrally (configurable).
- Outputs wrapped FASTA (60 chars/line).
- Very fast: uses os.urandom and byte-to-base mapping (no numpy needed).

Usage:
  python pad_to_10kb.py input.fasta output.fasta
  # optional:
  python pad_to_10kb.py input.fasta output.fasta --length 10000 --trim center --seed 42
"""

import argparse
import os
import sys

BASES = b"ACGT"

def rand_bases(n: int, seed: int = None) -> str:
    """Generate n random bases using os.urandom mapped modulo 4 (very fast)."""
    if n <= 0:
        return ""
    # os.urandom for speed; map bytes to 0..3 using modulo 4
    bs = os.urandom(n)
    # translate bytes to bases by modulo mapping
    # build a small lookup table once per call (fast enough)
    out = bytearray(n)
    for i, b in enumerate(bs):
        out[i] = BASES[b & 0x03]  # b % 4, "& 3" is faster
    return out.decode("ascii")

def wrap(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i:i+width]

def process_record(header: str, seq: str, target_len: int, trim_mode: str) -> str:
    L = len(seq)
    if L == target_len:
        # no change
        return seq
    if L > target_len:
        # trim
        if trim_mode == "center":
            start = (L - target_len) // 2
            return seq[start:start + target_len]
        elif trim_mode == "left":
            return seq[:target_len]
        elif trim_mode == "right":
            return seq[-target_len:]
        else:
            raise ValueError(f"Unknown trim_mode: {trim_mode}")
    # pad
    deficit = target_len - L
    left = deficit // 2
    right = deficit - left
    left_pad = rand_bases(left)
    right_pad = rand_bases(right)
    return left_pad + seq + right_pad

def run(in_fa: str, out_fa: str, target_len: int, trim_mode: str):
    with open(in_fa, "r") as fh, open(out_fa, "w") as fo:
        header = None
        buf = []
        for line in fh:
            if line.startswith(">"):
                # flush previous
                if header is not None:
                    seq = "".join(buf).replace(" ", "").replace("\n", "").replace("\r", "").upper()
                    seq = process_record(header, seq, target_len, trim_mode)
                    fo.write(header)
                    for chunk in wrap(seq, 60):
                        fo.write(chunk + "\n")
                header = line if line.endswith("\n") else (line + "\n")
                buf = []
            else:
                buf.append(line.strip())
        # flush last
        if header is not None:
            seq = "".join(buf).replace(" ", "").replace("\n", "").replace("\r", "").upper()
            seq = process_record(header, seq, target_len, trim_mode)
            fo.write(header)
            for chunk in wrap(seq, 60):
                fo.write(chunk + "\n")

def main():
    ap = argparse.ArgumentParser(description="Pad/truncate FASTA sequences to exactly N bp with random flanks.")
    ap.add_argument("input", help="Input FASTA")
    ap.add_argument("output", help="Output FASTA")
    ap.add_argument("--length", type=int, default=10000, help="Target length (default: 10000)")
    ap.add_argument("--trim", choices=["center", "left", "right"], default="center",
                    help="When longer than target, how to trim (default: center)")
    # Note: os.urandom is used; --seed kept for interface compatibility (no effect with os.urandom)
    ap.add_argument("--seed", type=int, default=None, help="(Ignored) Random seed placeholder")
    args = ap.parse_args()

    if args.length <= 0:
        sys.stderr.write("ERROR: --length must be positive\n")
        sys.exit(2)

    run(args.input, args.output, args.length, args.trim)

if __name__ == "__main__":
    main()
