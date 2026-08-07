#!/usr/bin/env python3
"""
make_genomic_context_benchmark_local.py

Replace random padding in the metapointfinder amino-acid point-mutation benchmark
with real genomic flanking sequence selected by LOCAL TBLASTX against a small
complete-genome BLAST database. Low-complexity masking is disabled by default
(-seg no / -dust no where applicable, -soft_masking false) and tblastx uses
-xdrop_ungap 900 by default to reduce split HSPs caused by ungapped extension cutoff.

This is intended as a faster, reproducible replacement for remote NCBI BLAST:
  1. Build a small local nucleotide BLAST database from complete RefSeq genomes.
  2. Use local tblastx to select a homologous genomic locus for each WT/RES pair.
  3. Insert both WT and RES alleles into the same genomic background.
  4. Output fixed-length pseudogenomes suitable for the existing wgsim step.

Example:
  python make_genomic_context_benchmark_local.py \
      benchmark/Eecoli.fasta \
      benchmark/Eecoli.6000.fasta \
      --db context_genomes/context_genomes.fna \
      --genomes context_genomes/context_genomes.fna \
      --length 6000 \
      --metadata benchmark/local_context_metadata.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BAD_TITLE_WORDS = [
    "plasmid",
    "contig",
    "scaffold",
    "whole genome shotgun",
    " wgs ",
    "draft genome",
    "partial sequence",
    "metagenome",
    "uncultured",
    "mag",
    "sag",
]

GOOD_TITLE_WORDS = [
    "complete genome",
    "complete sequence",
    "chromosome",
]


@dataclass
class LocalBlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    stitle: str

    @property
    def subject_start(self) -> int:
        return min(self.sstart, self.send)

    @property
    def subject_end(self) -> int:
        return max(self.sstart, self.send)

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

    @property
    def qcov(self) -> float:
        if self.qlen <= 0:
            return 0.0
        return (abs(self.qend - self.qstart) + 1) / self.qlen


@dataclass
class ContextInfo:
    group_key: str
    query_id: str
    subject_id: str
    subject_title: str
    pident: float
    qcov: float
    evalue: float
    bitscore: float
    strand: str
    hit_start: int
    hit_end: int
    subject_length: int
    fetch_start: int
    fetch_end: int
    insert_start_0based: int
    insert_end_0based: int
    allele_length: int
    final_length: int
    blast_program: str
    source: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use local tblastx hits against a small complete-genome BLAST database "
            "to place WT/RES benchmark alleles in the same real genomic context."
        )
    )
    parser.add_argument("input_fasta", help="Nucleotide FASTA produced after backtranseq")
    parser.add_argument("output_fasta", help="Output fixed-length genomic-context FASTA")
    parser.add_argument(
        "--db",
        required=True,
        help="Local nucleotide BLAST database prefix, usually context_genomes/context_genomes.fna",
    )
    parser.add_argument(
        "--genomes",
        required=True,
        help="Combined genome FASTA used to build the BLAST database; used to extract flanks",
    )
    parser.add_argument("--length", type=int, default=6000, help="Final pseudogenome length [6000]")
    parser.add_argument(
        "--blast-program",
        choices=["tblastx", "blastn"],
        default="tblastx",
        help="Local BLAST program for background selection [tblastx]",
    )
    parser.add_argument("--min-pident", type=float, default=50.0, help="Minimum percent identity [50]")
    parser.add_argument("--min-qcov", type=float, default=0.50, help="Minimum query coverage [0.50]")
    parser.add_argument("--evalue", default="1e-5", help="BLAST e-value threshold [1e-5]")
    parser.add_argument("--max-target-seqs", type=int, default=100, help="Max target sequences per query [100]")
    parser.add_argument("--num-threads", type=int, default=4, help="BLAST threads [4]")
    parser.add_argument(
        "--seg",
        choices=["yes", "no"],
        default="no",
        help="SEG low-complexity filtering for translated searches such as tblastx [no]",
    )
    parser.add_argument(
        "--dust",
        choices=["yes", "no"],
        default="no",
        help="DUST low-complexity filtering for nucleotide searches such as blastn [no]",
    )
    parser.add_argument(
        "--soft-masking",
        choices=["true", "false"],
        default="false",
        help="Use BLAST soft masking [false]",
    )
    parser.add_argument(
        "--xdrop-ungap",
        type=int,
        default=900,
        help="X-dropoff value for ungapped extensions, equivalent in spirit to old blastall -X [900]",
    )
    parser.add_argument(
        "--print-blast-command",
        action="store_true",
        help="Print each local BLAST command to stderr before running it",
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="Optional TSV file with selected context metadata",
    )
    parser.add_argument(
        "--cache",
        default=None,
        help="Optional JSON cache for selected contexts",
    )
    parser.add_argument(
        "--fail-on-missing",
        action="store_true",
        help=(
            "Stop with an error if a group has no suitable complete non-plasmid genome hit. "
            "By default, such groups are skipped and reported on the command line."
        ),
    )
    parser.add_argument(
        "--allow-bad-title",
        action="store_true",
        help="Do not reject hits based on plasmid/contig/scaffold/draft words in the title",
    )
    parser.add_argument(
        "--require-good-title",
        action="store_true",
        help="Require chromosome/complete-like wording in the BLAST subject title",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary query and BLAST output files for debugging",
    )
    return parser.parse_args()


def clean_id(record_id: str) -> str:
    return record_id.replace("|", "_").replace("=", "_")


def group_key_from_id(record_id: str) -> str:
    rid = clean_id(record_id)
    parts = rid.split("_")
    normalized = []
    replaced = False
    for p in parts:
        if p in {"WT", "RES"} and not replaced:
            normalized.append("ALLELE")
            replaced = True
        else:
            normalized.append(p)
    if replaced:
        return "_".join(normalized)
    rid = re.sub(r"_WT(_|$)", r"_ALLELE\1", rid)
    rid = re.sub(r"_RES(_|$)", r"_ALLELE\1", rid)
    return rid


def is_wt_record(record_id: str) -> bool:
    return "_WT" in clean_id(record_id).split()


def is_wt_id(record_id: str) -> bool:
    return "_WT" in clean_id(record_id)


def is_res_id(record_id: str) -> bool:
    return "_RES" in clean_id(record_id)


def load_records(path: str) -> List[SeqRecord]:
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    for r in records:
        r.id = clean_id(r.id)
        r.name = r.id
        r.description = r.id
        r.seq = Seq(str(r.seq).upper().replace("U", "T"))
    return records


def group_records(records: Iterable[SeqRecord]) -> Dict[str, List[SeqRecord]]:
    groups: Dict[str, List[SeqRecord]] = {}
    for record in records:
        groups.setdefault(group_key_from_id(record.id), []).append(record)
    return groups


def choose_query_record(group: List[SeqRecord]) -> SeqRecord:
    wt = [r for r in group if is_wt_id(r.id)]
    if wt:
        return sorted(wt, key=lambda r: r.id)[0]
    return sorted(group, key=lambda r: r.id)[0]


def title_is_bad(title: str) -> bool:
    t = f" {title.lower()} "
    return any(word in t for word in BAD_TITLE_WORDS)


def title_is_good(title: str) -> bool:
    t = title.lower()
    return any(word in t for word in GOOD_TITLE_WORDS)


def normalize_subject_id(value: str) -> str:
    """Return a set of practical BLAST/FASTA accession aliases as a single preferred string."""
    v = value.strip()
    # BLAST subject IDs often look like ref|NZ_CP007539.2| whereas FASTA IDs
    # may be NZ_CP007539.2, lcl|NZ_CP007539.2, or similar.
    if "|" in v:
        parts = [x for x in v.split("|") if x]
        for part in parts:
            if re.match(r"^[A-Z]{1,3}_[A-Z]*[0-9]+\.[0-9]+$", part) or re.match(r"^[A-Z]{1,4}[0-9]+\.[0-9]+$", part):
                return part
        if parts:
            return parts[-1]
    return v


def subject_aliases(value: str) -> List[str]:
    aliases = []
    v = value.strip()
    aliases.append(v)
    aliases.append(normalize_subject_id(v))
    if "|" in v:
        aliases.extend([x for x in v.split("|") if x])
    # Remove duplicate aliases while preserving order.
    seen = set()
    out = []
    for a in aliases:
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return out


def read_genomes(genome_fasta: str) -> Dict[str, SeqRecord]:
    records: Dict[str, SeqRecord] = {}
    for rec in SeqIO.parse(genome_fasta, "fasta"):
        # Store the original ID plus common accession aliases so BLAST IDs and
        # FASTA IDs can be matched robustly.
        for alias in subject_aliases(rec.id):
            records.setdefault(alias, rec)
        if rec.description:
            for alias in subject_aliases(rec.description.split()[0]):
                records.setdefault(alias, rec)
    if not records:
        raise ValueError(f"No records found in genome FASTA: {genome_fasta}")
    return records


def run_local_blast(record: SeqRecord, args: argparse.Namespace, tmpdir: Path) -> List[LocalBlastHit]:
    query = tmpdir / f"query_{record.id}.fasta"
    out = tmpdir / f"query_{record.id}.blast.tsv"
    SeqIO.write([record], query, "fasta")

    outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle"
    cmd = [
        args.blast_program,
        "-query", str(query),
        "-db", args.db,
        "-out", str(out),
        "-outfmt", outfmt,
        "-evalue", str(args.evalue),
        "-max_target_seqs", str(args.max_target_seqs),
        "-num_threads", str(args.num_threads),
        "-soft_masking", args.soft_masking,
        "-xdrop_ungap", str(args.xdrop_ungap),
    ]
    if args.blast_program == "tblastx":
        # tblastx is an ungapped translated search in BLAST+.
        # Use SEG, not DUST, for translated low-complexity filtering.
        cmd.extend(["-seg", args.seg])
    elif args.blast_program == "blastn":
        # DUST applies to nucleotide BLAST searches.
        cmd.extend(["-dust", args.dust])

    if args.print_blast_command:
        print("BLASTCMD\t" + " ".join(cmd), file=sys.stderr)

    try:
        subprocess.run(cmd, check=True, text=True, capture_output=True)
    except FileNotFoundError as e:
        raise RuntimeError(f"Could not find local BLAST executable '{args.blast_program}' in PATH") from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Local BLAST failed for {record.id}\nCommand: {' '.join(cmd)}\nSTDERR:\n{e.stderr}"
        ) from e

    hits: List[LocalBlastHit] = []
    if not out.exists() or out.stat().st_size == 0:
        return hits
    with out.open() as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 13:
                continue
            hits.append(
                LocalBlastHit(
                    qseqid=row[0],
                    sseqid=row[1],
                    pident=float(row[2]),
                    length=int(row[3]),
                    qstart=int(row[4]),
                    qend=int(row[5]),
                    sstart=int(row[6]),
                    send=int(row[7]),
                    evalue=float(row[8]),
                    bitscore=float(row[9]),
                    qlen=int(row[10]),
                    slen=int(row[11]),
                    stitle=row[12],
                )
            )
    return hits


def passes_hit_filters(hit: LocalBlastHit, args: argparse.Namespace) -> Tuple[bool, str]:
    if hit.pident < args.min_pident:
        return False, f"pident {hit.pident:.2f} < {args.min_pident:.2f}"
    if hit.qcov < args.min_qcov:
        return False, f"qcov {hit.qcov:.3f} < {args.min_qcov:.3f}"
    if not args.allow_bad_title and title_is_bad(hit.stitle):
        return False, "bad subject title"
    if args.require_good_title and not title_is_good(hit.stitle):
        return False, "subject title lacks complete/chromosome wording"
    if hit.subject_start < 1 or hit.subject_end > hit.slen:
        return False, "subject coordinates outside subject sequence"
    return True, "OK"


def select_best_hit(hits: List[LocalBlastHit], args: argparse.Namespace) -> Optional[LocalBlastHit]:
    accepted: List[LocalBlastHit] = []
    for hit in hits:
        ok, _ = passes_hit_filters(hit, args)
        if ok:
            accepted.append(hit)
    if not accepted:
        return None
    return sorted(accepted, key=lambda h: (h.bitscore, h.pident, h.qcov, -h.evalue), reverse=True)[0]


def fit_context_window(hit: LocalBlastHit, allele_len: int, final_len: int) -> Tuple[int, int, int, int]:
    if allele_len > final_len:
        raise ValueError(f"Allele length {allele_len} exceeds final length {final_len}")

    hit_start = hit.subject_start
    hit_end = hit.subject_end
    hit_center = (hit_start + hit_end) // 2

    fetch_start = hit_center - final_len // 2
    fetch_end = fetch_start + final_len - 1

    if fetch_start < 1:
        fetch_start = 1
        fetch_end = final_len
    if fetch_end > hit.slen:
        fetch_end = hit.slen
        fetch_start = hit.slen - final_len + 1

    if fetch_start < 1 or fetch_end > hit.slen or (fetch_end - fetch_start + 1) != final_len:
        raise ValueError(
            f"Cannot fit {final_len} bp context around hit {hit.sseqid}:{hit_start}-{hit_end}; subject length={hit.slen}"
        )

    insert_center = hit_center - fetch_start
    insert_start = insert_center - allele_len // 2
    insert_end = insert_start + allele_len

    if insert_start < 0:
        insert_start = 0
        insert_end = allele_len
    if insert_end > final_len:
        insert_end = final_len
        insert_start = final_len - allele_len

    return fetch_start, fetch_end, insert_start, insert_end


def extract_window(genomes: Dict[str, SeqRecord], subject_id: str, start: int, end: int) -> SeqRecord:
    rec = None
    for alias in subject_aliases(subject_id):
        if alias in genomes:
            rec = genomes[alias]
            break
    if rec is None:
        raise KeyError(
            f"Subject {subject_id} not found in genome FASTA. Tried aliases: {', '.join(subject_aliases(subject_id))}"
        )
    seq = rec.seq[start - 1:end]
    return SeqRecord(seq, id=f"{rec.id}:{start}-{end}", description=rec.description)


def make_context_sequence(background: SeqRecord, allele: SeqRecord, insert_start: int, insert_end: int, strand: str) -> SeqRecord:
    bg_seq = background.seq
    allele_seq = allele.seq
    if strand == "-":
        allele_seq = allele_seq.reverse_complement()
    new_seq = bg_seq[:insert_start] + allele_seq + bg_seq[insert_end:]
    return SeqRecord(Seq(str(new_seq).upper()), id=allele.id, name=allele.id, description=allele.id)


def load_cache(path: Optional[str]) -> Dict[str, dict]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    return json.loads(p.read_text())


def save_cache(path: Optional[str], cache: Dict[str, dict]) -> None:
    if not path:
        return
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(cache, indent=2, sort_keys=True))


def cache_key(group_key: str, args: argparse.Namespace) -> str:
    return (
        f"local:{args.blast_program}:{args.db}:{args.length}:{args.min_pident}:{args.min_qcov}:"
        f"seg={args.seg}:dust={args.dust}:soft={args.soft_masking}:xdrop={args.xdrop_ungap}:"
        f"{group_key}"
    )


def context_from_cache(group_key: str, cache: Dict[str, dict], args: argparse.Namespace) -> Optional[ContextInfo]:
    value = cache.get(cache_key(group_key, args))
    if not value:
        return None
    return ContextInfo(**value)


def put_context_cache(group_key: str, info: ContextInfo, cache: Dict[str, dict], args: argparse.Namespace) -> None:
    cache[cache_key(group_key, args)] = asdict(info)


def metadata_header() -> List[str]:
    return list(ContextInfo.__dataclass_fields__.keys())


def write_metadata(path: Optional[str], rows: List[ContextInfo]) -> None:
    if not path:
        return
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_header(), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def find_context_for_group(group_key: str, group: List[SeqRecord], args: argparse.Namespace, tmpdir: Path, cache: Dict[str, dict]) -> ContextInfo:
    cached = context_from_cache(group_key, cache, args)
    if cached:
        return cached

    query = choose_query_record(group)
    hits = run_local_blast(query, args, tmpdir)
    if not hits:
        raise RuntimeError("no local BLAST hits returned")
    hit = select_best_hit(hits, args)
    if hit is None:
        reasons = []
        for h in hits[:10]:
            _, reason = passes_hit_filters(h, args)
            reasons.append(f"{h.sseqid}:{reason}")
        raise RuntimeError("no suitable hit after filtering; first reasons: " + "; ".join(reasons))

    fetch_start, fetch_end, insert_start, insert_end = fit_context_window(hit, len(query.seq), args.length)
    info = ContextInfo(
        group_key=group_key,
        query_id=query.id,
        subject_id=hit.sseqid,
        subject_title=hit.stitle,
        pident=hit.pident,
        qcov=hit.qcov,
        evalue=hit.evalue,
        bitscore=hit.bitscore,
        strand=hit.strand,
        hit_start=hit.subject_start,
        hit_end=hit.subject_end,
        subject_length=hit.slen,
        fetch_start=fetch_start,
        fetch_end=fetch_end,
        insert_start_0based=insert_start,
        insert_end_0based=insert_end,
        allele_length=len(query.seq),
        final_length=args.length,
        blast_program=args.blast_program,
        source="local_blast",
    )
    put_context_cache(group_key, info, cache, args)
    return info


def process(args: argparse.Namespace) -> int:
    if shutil.which(args.blast_program) is None:
        raise RuntimeError(f"Could not find {args.blast_program} in PATH")
    records = load_records(args.input_fasta)
    groups = group_records(records)
    genomes = read_genomes(args.genomes)
    cache = load_cache(args.cache)

    output_records: List[SeqRecord] = []
    metadata: List[ContextInfo] = []
    skipped: List[Tuple[str, str, str]] = []

    tmp_parent = None if args.keep_temp else tempfile.TemporaryDirectory(prefix="mpf_local_context_")
    tmpdir = Path(tmp_parent.name if tmp_parent else tempfile.mkdtemp(prefix="mpf_local_context_keep_"))

    try:
        for group_key in sorted(groups):
            group = groups[group_key]
            query = choose_query_record(group)
            try:
                info = find_context_for_group(group_key, group, args, tmpdir, cache)
                background = extract_window(genomes, info.subject_id, info.fetch_start, info.fetch_end)
                for allele in sorted(group, key=lambda r: r.id):
                    out = make_context_sequence(
                        background,
                        allele,
                        info.insert_start_0based,
                        info.insert_end_0based,
                        info.strand,
                    )
                    if len(out.seq) != args.length:
                        raise RuntimeError(f"Output length for {out.id} is {len(out.seq)}, expected {args.length}")
                    output_records.append(out)
                metadata.append(info)
                print(
                    f"OK\tgroup={group_key}\tquery={query.id}\tsubject={info.subject_id}\t"
                    f"pident={info.pident:.1f}\tqcov={info.qcov:.3f}\tstrand={info.strand}",
                    file=sys.stderr,
                )
            except Exception as e:
                msg = str(e)
                if args.fail_on_missing:
                    raise
                skipped.append((group_key, query.id, msg))
                print(f"SKIP\tgroup={group_key}\tquery={query.id}\treason={msg}", file=sys.stderr)

        if not output_records:
            raise RuntimeError("No output records were created; all groups were skipped")
        Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(output_records, args.output_fasta, "fasta")
        write_metadata(args.metadata, metadata)
        save_cache(args.cache, cache)

        print(
            f"Wrote {len(output_records)} records from {len(metadata)} context groups to {args.output_fasta}. "
            f"Skipped {len(skipped)} of {len(groups)} groups. BLAST program: {args.blast_program}",
            file=sys.stderr,
        )
        if skipped:
            print("Skipped groups:", file=sys.stderr)
            for g, q, reason in skipped:
                print(f"SKIPPED\t{g}\t{q}\t{reason}", file=sys.stderr)
        return 0
    finally:
        if tmp_parent is not None:
            tmp_parent.cleanup()
        else:
            print(f"Temporary files kept in: {tmpdir}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    return process(args)


if __name__ == "__main__":
    raise SystemExit(main())
