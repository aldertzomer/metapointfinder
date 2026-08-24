#!/usr/bin/env python3
"""Score DNA mutations from a tab-separated MetaPointFinder input file.

Python translation of ``dna_score_mutations.R``.  The reference is aligned
globally to a local region of the target, using the same match, mismatch and
affine-gap scores as the R implementation.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path


OUTPUT_FILE = "updated_table_with_scores_and_mutations.tsv"
TOKEN_RE = re.compile(r"^([ACGT]+)([0-9]+)([ACGT-]+)$")
NEG_INF = -math.inf


@dataclass(frozen=True)
class Mutation:
    ref: str
    position: int
    alt: str
    raw: str


def _global_local_align(reference: str, target: str) -> tuple[str, str]:
    """Globally align reference to a local target region with affine gaps."""
    m, n = len(reference), len(target)
    # M: paired bases; X: reference base/target gap; Y: reference gap/target base.
    matrices = [[[NEG_INF] * (n + 1) for _ in range(m + 1)] for _ in range(3)]
    traces = [[[-1] * (n + 1) for _ in range(m + 1)] for _ in range(3)]
    match, mismatch, gap_open, gap_extend = 2.0, -5.0, -10.0, -0.5

    matrices[0][0][0] = 0.0
    # A target prefix is free (local with respect to the subject).
    for j in range(1, n + 1):
        matrices[0][0][j] = 0.0
        traces[0][0][j] = 0
    for i in range(1, m + 1):
        matrices[1][i][0] = gap_open + (i - 1) * gap_extend
        traces[1][i][0] = 0 if i == 1 else 1

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            candidates = [matrices[s][i - 1][j - 1] for s in range(3)]
            previous = max(range(3), key=candidates.__getitem__)
            score = match if reference[i - 1] == target[j - 1] else mismatch
            matrices[0][i][j] = candidates[previous] + score
            traces[0][i][j] = previous

            candidates = [
                matrices[0][i - 1][j] + gap_open,
                matrices[1][i - 1][j] + gap_extend,
                matrices[2][i - 1][j] + gap_open,
            ]
            previous = max(range(3), key=candidates.__getitem__)
            matrices[1][i][j] = candidates[previous]
            traces[1][i][j] = previous

            candidates = [
                matrices[0][i][j - 1] + gap_open,
                matrices[1][i][j - 1] + gap_open,
                matrices[2][i][j - 1] + gap_extend,
            ]
            previous = max(range(3), key=candidates.__getitem__)
            matrices[2][i][j] = candidates[previous]
            traces[2][i][j] = previous

    state, end_j = max(
        ((s, j) for s in range(3) for j in range(n + 1)),
        key=lambda item: matrices[item[0]][m][item[1]],
    )
    i, j = m, end_j
    aligned_ref: list[str] = []
    aligned_target: list[str] = []
    while i > 0:
        previous = traces[state][i][j]
        if state == 0:
            aligned_ref.append(reference[i - 1])
            aligned_target.append(target[j - 1])
            i -= 1
            j -= 1
        elif state == 1:
            aligned_ref.append(reference[i - 1])
            aligned_target.append("-")
            i -= 1
        else:
            aligned_ref.append("-")
            aligned_target.append(target[j - 1])
            j -= 1
        state = previous
    return "".join(reversed(aligned_ref)), "".join(reversed(aligned_target))


def _unknown() -> dict[str, object]:
    return {
        "score": 0,
        "detected_mutations": "None",
        "wt_confirmed_positions": 0,
        "status": "Unknown",
    }


def calculate_mutation_score(
    reference: str | None, target: str | None, changes_str: str | None
) -> dict[str, object]:
    if not changes_str or changes_str.strip().upper() in {"", "NA", "NAN"}:
        return _unknown()

    mutations: list[Mutation] = []
    for token in changes_str.split(","):
        raw = token.strip().upper()
        match = TOKEN_RE.fullmatch(raw)
        if match:
            mutations.append(Mutation(match[1], int(match[2]), match[3], raw))
    if not mutations:
        return _unknown()

    reference = re.sub(r"\s", "", reference or "").upper()
    target = re.sub(r"[-\s]", "", target or "").upper()
    if not reference or not target:
        return _unknown()
    aligned_ref, aligned_target = _global_local_align(reference, target)

    ref_positions: list[int | None] = []
    ref_position = 0
    for base in aligned_ref:
        if base != "-":
            ref_position += 1
            ref_positions.append(ref_position)
        else:
            ref_positions.append(None)

    covered = [i for i, base in enumerate(aligned_target) if base != "-"]
    if not covered:
        return _unknown()
    coverage_start, coverage_end = min(covered), max(covered)
    score = 0
    detected: list[str] = []
    detected_positions: set[int] = set()
    wt_positions: set[int] = set()

    def add_wt(position: int) -> None:
        if position not in detected_positions:
            wt_positions.add(position)

    for mutation in mutations:
        pos, ref_seg, alt_seg = mutation.position, mutation.ref, mutation.alt
        ref_len, alt_len = len(ref_seg), len(alt_seg)
        if pos < 1 or pos + ref_len - 1 > ref_position:
            continue
        wanted = set(range(pos, pos + ref_len))
        ref_cols = [i for i, value in enumerate(ref_positions) if value in wanted]
        if (len(ref_cols) != ref_len or min(ref_cols) < coverage_start
                or max(ref_cols) > coverage_end):
            continue
        ref_block = "".join(aligned_ref[i] for i in ref_cols)
        read_block = "".join(aligned_target[i] for i in ref_cols)
        if ref_block != ref_seg:
            continue

        found = False
        if alt_seg == "-":
            found = all(aligned_target[i] == "-" for i in ref_cols)
            if not found and read_block == ref_seg and "N" not in read_block:
                add_wt(pos)
        elif ref_len == 1 and alt_len > 1 and alt_seg.startswith(ref_seg):
            col = ref_cols[0]
            base_ok = aligned_target[col] == ref_seg
            needed = alt_seg[1:]
            seen: list[str] = []
            started = False
            next_col = col + 1
            while next_col < len(aligned_ref) and len(seen) < len(needed):
                if aligned_ref[next_col] == "-":
                    base = aligned_target[next_col]
                    if base not in {"-", "N"}:
                        seen.append(base)
                        started = True
                elif started:
                    break
                next_col += 1
            found = base_ok and "".join(seen) == needed
            if not found and base_ok:
                add_wt(pos)
        elif ref_len == alt_len and alt_seg != ref_seg:
            found = "-" not in read_block and "N" not in read_block and read_block == alt_seg
            if not found and "-" not in read_block and "N" not in read_block and read_block == ref_seg:
                add_wt(pos)
        elif "-" not in read_block and "N" not in read_block and read_block == ref_seg:
            add_wt(pos)

        if found:
            score += 1
            detected.append(ref_seg + str(pos) + "-" if alt_seg == "-" else mutation.raw)
            detected_positions.add(pos)
            wt_positions.discard(pos)

    return {
        "score": score,
        "detected_mutations": ",".join(detected) if detected else "None",
        "wt_confirmed_positions": len(wt_positions),
        "status": "Resistant" if score else ("Wildtype" if wt_positions else "Unknown"),
    }


def _score_row(row: dict[str, str]) -> dict[str, object]:
    return calculate_mutation_score(row["reference"], row["target"], row["changes_str"])


def process_mutation_data(file_path: str | Path, threads: int = 1) -> list[dict[str, str]]:
    with Path(file_path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        required = {"class", "gene", "read", "reference", "target", "changes_str"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    if threads > 1 and len(rows) > 1:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            results = list(executor.map(_score_row, rows))
    else:
        results = [_score_row(row) for row in rows]

    output_columns = [
        "MutationScore", "DetectedMutations", "WTConfirmedPositions", "Status"
    ]
    for row, result in zip(rows, results):
        row.update(dict(zip(output_columns, result.values())))
    with Path(OUTPUT_FILE).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames + output_columns, delimiter="\t",
            quoting=csv.QUOTE_NONE, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", help="Input TSV file")
    parser.add_argument("threads", nargs="?", type=int, default=1)
    args = parser.parse_args()
    rows = process_mutation_data(args.input_file, max(1, args.threads))
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
