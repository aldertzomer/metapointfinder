#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
metapointfinder.py
Reimplementation of the Bash metapointfinder pipeline in Python

Usage:
  python metapointfinder.py --input <file.fastq[.gz]> --db <databasefolder> --output <outputfolder> --id 85 --threads 4 [--force]

Requirements:
  - Programs listed in ./dependencies (diamond, kma, wget, R)
  - R packages: Biostrings, pwalign, parallel
  - DIAMOND 2.0.15

This script expects to live next to:
  - dependencies
  - prot_score_mutations.R
  - dna_score_mutations.R
"""

import os
import sys
import re
import shutil
import gzip
import subprocess
from pathlib import Path

#
# -------- Add helper so class summaries include DB classes even if 0/0/0
#
ALL_CLASSES = set()

def collect_all_classes(db_dir: Path):
    classes = set()
    # Protein mutation TSV: class typically at column 6
    prot_mut = Path(db_dir) / "AMRProt-mutation.tsv"
    if prot_mut.exists():
        with open(prot_mut, "r", encoding="utf-8") as f:
            for ln in f:
                if not ln.strip():
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 6 and parts[5] and parts[5] != "class":
                    classes.add(parts[5])
    # DNA mutation TSV: class typically at column 5
    dna_mut = Path(db_dir) / "AMR_DNA-mutation.tsv"
    if dna_mut.exists():
        with open(dna_mut, "r", encoding="utf-8") as f:
            for ln in f:
                if not ln.strip():
                    continue
                if "mutation_position" in ln or ("class" in ln and "position" in ln):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 5 and parts[4] and parts[4] != "class":
                    classes.add(parts[4])
    return classes

# ------------------------
# Small helpers
# ------------------------

def sh(cmd, cwd=None, shell=True, stdout=None, stderr=None):
    """Run a shell command, raise on failure."""
    return subprocess.run(cmd, cwd=cwd, shell=shell, check=True,
                          stdout=stdout or subprocess.PIPE,
                          stderr=stderr or subprocess.PIPE)

def exists(p):
    return Path(p).exists()

def read_text(p):
    return Path(p).read_text(encoding="utf-8")

def write_text(p, s):
    Path(p).write_text(s, encoding="utf-8")

def ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)

def which_or_die(prog):
    if shutil.which(prog) is None:
        print(f"Required dependency {prog} not installed")
        sys.exit(1)

def fasta_to_tsv_two_cols(fa_in: Path, tsv_out: Path):
    """
    Convert FASTA into two-column TSV: accession<TAB>sequence
    Accession is the header up to first whitespace
    """
    acc = None
    seq_parts = []
    with open(fa_in, "r", encoding="utf-8") as fin, open(tsv_out, "w", encoding="utf-8") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    fout.write(f"{acc}\t{''.join(seq_parts)}\n")
                acc = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if acc is not None:
            fout.write(f"{acc}\t{''.join(seq_parts)}\n")

# ------------------------
# Mutation combiners for database - future versions will include predicted similars based on structural similarities (e.g. S83I looks a lot like S83L in GyrA and should be included as Resistant marker)
# ------------------------

def parse_amrprot_mutation_underscore_combined(tsv_in):
    """
    Build dict: (class, accession) -> sorted, comma-joined changes_str (low->high).
      - Use the 'position' column (numeric) from AMRProt-mutation.tsv (field 2,3,4,6)
      - Keep raw token until the END, then map 'del'/'STOP' -> '-'
      - Do NOT strip '-' from the right side (so '.107-' is preserved)
      - Sort by position ascending
    """
    tmp = {}  # (class, accession) -> list[(pos_int, mutation_raw_string)]
    with open(tsv_in, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            # Expect the original columns; we mimic cut -f 2,3,4,6 path
            # Guard against header
            if "class" in line and "mutation" in line:
                continue
            if len(parts) < 6:
                # Best effort fallback when columns not in expected width
                continue
            accession = parts[1]
            position  = parts[2]
            mutation_original = parts[3]
            amr_class = parts[5]

            # Keep raw token; don't map del/STOP yet
            if "_" in mutation_original:
                mutation_cleaned = mutation_original.split("_", 1)[1]
            else:
                mutation_cleaned = mutation_original

            mpos = re.search(r"[0-9]+", position)
            if not mpos:
                continue
            pos_int = int(mpos.group(0))

            # Split around the digits; keep any '-' in baseright
            pieces = re.split(r"([0-9]+)", mutation_cleaned, maxsplit=1)
            if len(pieces) < 3:
                continue
            baseleft  = pieces[0]
            baseright = pieces[2]

            mutation_raw = f"{baseleft}{pos_int}{baseright}"
            tmp.setdefault((amr_class, accession), []).append((pos_int, mutation_raw))

    combined = {}
    for key, items in tmp.items():
        items.sort(key=lambda x: x[0])  # low -> high
        joined = ",".join(m for _, m in items)
        joined = joined.replace("del", "-").replace("STOP", "-")
        combined[key] = joined
    return combined

def build_AMR_DNA_combined(db_dir):
    """
    - Always use the 'position in provided sequence' column (parts[3]) as the coordinate
      when available; otherwise, fall back to digits from parts[1].
    - Strip trailing '-' from the left token segment (handles promoter-style negatives like G-48T).
    """
    from pathlib import Path
    dna_mut       = Path(db_dir) / "AMR_DNA-mutation.tsv"
    dna_mut_us    = Path(db_dir) / "AMR_DNA-mutation_underscore.tsv"
    dna_fa        = Path(db_dir) / "AMR_DNA.fa"
    dna_fa_us     = Path(db_dir) / "AMR_DNA_underscore.fa"
    dna_fa_tsv    = Path(db_dir) / "AMR_DNA_underscore_fa.tsv"
    combined_out  = Path(db_dir) / "AMR_DNA_underscore_combined.tsv"

    if combined_out.exists():
        return

    # underscore replacements
    write_text(dna_mut_us, read_text(dna_mut).replace("/", "_"))
    write_text(dna_fa_us,  read_text(dna_fa).replace("/", "_"))

    # FASTA -> TSV
    fasta_to_tsv_two_cols(dna_fa_us, dna_fa_tsv)
    acc_to_seq = {}
    with open(dna_fa_tsv, "r", encoding="utf-8") as f:
        for line in f:
            acc, seq = line.rstrip("\n").split("\t", 1)
            acc_to_seq[acc] = seq

    # Read mutation rows
    rows = []  # (acc, pos_used:int, mut_token, class)
    with open(dna_mut_us, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            if "mutation_position" in line or ("class" in line and "position" in line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            acc = parts[0]
            pos_str = parts[1]
            mut = parts[2]
            cls = parts[4]

            # ALWAYS prefer the override column (position in provided sequence)
            pos_used = None
            if len(parts) > 3 and parts[3].strip().isdigit():
                pos_used = int(parts[3].strip())
            else:
                mpos = re.search(r"[0-9]+", pos_str)
                if not mpos:
                    continue
                pos_used = int(mpos.group(0))

            rows.append((acc, pos_used, mut, cls))

    classes    = sorted(set(r[3] for r in rows))
    accessions = sorted(set(r[0] for r in rows))

    by_cls_acc = {}
    for cls in classes:
        for acc in accessions:
            muts = []
            for (acc_i, pos_used, mut, cls_i) in rows:
                if acc_i != acc or cls_i != cls:
                    continue

                mutation_cleaned = mut.split("_", 1)[1] if "_" in mut else mut

                # Split token and strip trailing '-' from left part (promoter negative marker),
                # but keep '-' on the right (true deletion marker).
                pieces = re.split(r"([0-9]+)", mutation_cleaned, maxsplit=1)
                if len(pieces) < 3:
                    continue
                baseleft  = pieces[0].rstrip('-')
                baseright = pieces[2]

                mutation_raw = f"{baseleft}{pos_used}{baseright}"
                muts.append((pos_used, mutation_raw))

            if muts:
                muts.sort(key=lambda x: x[0])
                changes_str = ",".join(m for _, m in muts)
                changes_str = changes_str.replace("del", "-").replace("STOP", "-")
                seq = acc_to_seq.get(acc, "")
                by_cls_acc[(cls, acc)] = (seq, changes_str)

    # Write sorted by (class, accession)
    with open(combined_out, "w", encoding="utf-8") as out:
        for (cls, acc) in sorted(by_cls_acc.keys(), key=lambda t: (t[0], t[1])):
            seq, changes = by_cls_acc[(cls, acc)]
            out.write(f"{cls}\t{acc}\t{seq}\t{changes}\n")


# ------------------------
# Summaries
# ------------------------

def tsv_read(path):
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return header, rows

def write_class_and_gene_summaries(sample_prefix, out_dir: Path, is_prot=True):
    """
    Replicates:
      class summary -> <sample>.class.[prot|dna].summary.txt
      gene  summary -> <sample>.gene.[prot|dna].summary.txt
    from *.[prot|dna].updated_table_with_scores_and_mutations.tsv
    """
    kind = "prot" if is_prot else "dna"
    infile = out_dir / f"{sample_prefix}.{kind}.updated_table_with_scores_and_mutations.tsv"
    if not infile.exists():
        return
    header, rows = tsv_read(infile)

    # columns: class, gene, Status
    row_classes = sorted(set(r["class"] for r in rows if r.get("class") and r["class"] != "class"))
    # include DB classes as well (print zero rows too)
    classes = sorted(set(row_classes) | set(ALL_CLASSES))
    genes   = sorted(set(r["gene"]  for r in rows if r.get("gene")  and r["gene"]  != "gene" ))

    # class summary
    class_out = out_dir / f"{sample_prefix}.class.{kind}.summary.txt"
    with open(class_out, "w", encoding="utf-8") as f:
        f.write("class\tWT\tR\tUNKNOWN\n")
        for c in classes:
            WT = sum(1 for r in rows if r["class"] == c and r["Status"] == "Wildtype")
            R  = sum(1 for r in rows if r["class"] == c and r["Status"] == "Resistant")
            NA = sum(1 for r in rows if r["class"] == c and r["Status"] == "Unknown")
            f.write(f"{c}\t{WT}\t{R}\t{NA}\n")

    # gene summary
    gene_out = out_dir / f"{sample_prefix}.gene.{kind}.summary.txt"
    with open(gene_out, "w", encoding="utf-8") as f:
        f.write("gene\tWT\tR\tUNKNOWN\n")
        for c in classes:
            for g in genes:
                WT = sum(1 for r in rows if r["class"] == c and r["gene"] == g and r["Status"] == "Wildtype")
                R  = sum(1 for r in rows if r["class"] == c and r["gene"] == g and r["Status"] == "Resistant")
                NA = sum(1 for r in rows if r["class"] == c and r["gene"] == g and r["Status"] == "Unknown")
                if WT == 0 and R == 0 and NA == 0:
                    continue
                f.write(f"{c}.{g}\t{WT}\t{R}\t{NA}\n")

# ------------------------
# Main
# ------------------------

def main():
    print("This is metapointfinder v0.3")
    print("This tool finds substitions in translated reads matching to relevant proteins and mutations in reads matching to relevant genes")
    print("Relevant proteins and genes are obtained from the AMRFinder database")
    #
    # CLI flags shim: support --input/--db/--output [--force] while preserving positional form
    #
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--input", "-i")
    parser.add_argument("--db", "-d")
    parser.add_argument("--output", "-o")
    parser.add_argument("--identity", "-id")
    parser.add_argument("--threads", "-t")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("rest", nargs="*")
    args = parser.parse_args()

    if args.identity is None:
     args.identity = 85  # default identity threshold (1–100 scale)

    if args.threads is None:
     args.threads = 4  # default number of threads


    if args.input and args.db and args.output:
        input_path  = Path(args.input).resolve()
        database    = Path(args.db)
        output      = Path(args.output)
        identity    = args.identity
        threads     = args.threads
        force_flag  = "force" if args.force else None
    else:
        if len(args.rest) < 4:
            print("usage: metapointfinder.py --input file.fastq[.gz] --db databasefolder --output outputfolder --identity 85 --threads 4 [--force]")
#            print("   or: metapointfinder.py file.fastq[.gz] databasefolder outputfolder identity [force]")
            sys.exit(1)
#        input_path  = Path(args.rest[0]).resolve()
#        database    = Path(args.rest[1])
#        output      = Path(args.rest[2])
#        identity    = Path(args.rest[3])
#        force_flag  = "force" if (len(args.rest) > 3 and args.rest[3] == "force") else None

    script_path = Path(__file__).resolve()
    script_dir  = script_path.parent
    fastqfile = input_path.name

    # ------------------------
    # Dependencies (from ./dependencies)
    # ------------------------
    deps_file = script_dir / "dependencies"
    if not deps_file.exists():
        print("Missing 'dependencies' file next to the script.")
        sys.exit(1)

    # Check each dependency
    with open(deps_file, "r", encoding="utf-8") as f:
        for program in [ln.strip() for ln in f if ln.strip()]:
            if shutil.which(program) is None:
                print(f"Required dependency {program} not installed")
                sys.exit(1)

    # R packages
    try:
        sh('R -e "stopifnot(3 == length(find.package(c(\'parallel\',\'Biostrings\', \'pwalign\'))))" --slave')
    except subprocess.CalledProcessError:
        print("R package parallel and/or Biostrings and/or pwalign not installed")
        sys.exit(1)

    # Input checks
    if not input_path.exists():
        print(f"File {input_path} not found!")
        print("usage: metapointfinder.py file.fastq[.gz] databasefolder outputfolder [force]")
        sys.exit(1)

    # Database folder
    if database.exists():
        print("Database folder already exists")
        database = database.resolve()
    else:
        ensure_dir(database)
        database = database.resolve()
    if not database.exists():
        print("Error: database folder not existing or not created. Try again - exiting")
        sys.exit(1)

    # Output folder (respect 'force')
    if output.exists():
        if force_flag == "force":
            output = output.resolve()
        else:
            print("Error: output folder already exists.")
            sys.exit(1)
    else:
        ensure_dir(output)
        output = output.resolve()
    if not output.exists():
        print("Error: output folder not existing or not created. Try again - exiting")
        sys.exit(1)

    # Unzip/prepare sample FASTQ in output
    if fastqfile.endswith(".gz"):
        sample = fastqfile.replace(".fastq.gz", "")
        with gzip.open(input_path, "rb") as fin, open(output / f"{sample}.fastq", "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        sample = fastqfile.replace(".fastq", "")
        shutil.copyfile(input_path, output / f"{sample}.fastq")

    # ------------------------
    # Database prep
    # ------------------------
    os.chdir(database)
    if not exists("AMRProt.fa"):
        print("downloading protein mutation sequences")
        base = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
        sh(f"wget -nv -c {base}/AMRProt-mutation.tsv 1>> databaseprep.log 2>> databaseprep.error")
        sh(f"wget -nv -c {base}/AMRProt.fa 1>> databaseprep.log 2>> databaseprep.error")
    if not exists("AMRProt-mutation_underscore.tsv"):
        write_text("AMRProt-mutation_underscore.tsv", read_text("AMRProt-mutation.tsv").replace("/", "_"))

    # Combined protein changes
    if not exists("AMRProt-mutation_underscore_combined.tsv"):
        print("Prepping AMRProt mutation database")
        combined = parse_amrprot_mutation_underscore_combined("AMRProt-mutation_underscore.tsv")
        with open("AMRProt-mutation_underscore_combined.tsv", "w", encoding="utf-8") as out:
            for (cls, acc), changes in sorted(combined.items(), key=lambda t: (t[0][0], t[0][1])):
                out.write(f"{cls}\t{acc}\t{changes}\n")

    # Build AMRProt_mutation.fa (only headers containing "mutation")
    if not exists("AMRProt_mutation.fa"):
        lines = []
        with open("AMRProt.fa", "r", encoding="utf-8") as f:
            buf = f.read()
        # Split records
        recs = [r for r in buf.split(">") if r]
        for rec in recs:
            if "mutation" in rec:
                lines.append(">" + rec)
        write_text("AMRProt_mutation.fa", "".join(lines))

    # DIAMOND DB (only if missing)
    if not exists("AMRProt.dmnd"):
        print("Building DIAMOND database for AMRProt_mutation.fa")
        sh("diamond makedb --in AMRProt_mutation.fa -d AMRProt 1>> databaseprep.log 2>> databaseprep.error")

    if not exists("AMR_DNA-mutation.tsv"):
        print("downloading AMR_DNA-*.tsv")
        base = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
        dna_tsvs = [
            "AMR_DNA-Acinetobacter_baumannii.tsv",
            "AMR_DNA-Campylobacter.tsv",
            "AMR_DNA-Clostridioides_difficile.tsv",
            "AMR_DNA-Enterococcus_faecalis.tsv",
            "AMR_DNA-Enterococcus_faecium.tsv",
            "AMR_DNA-Escherichia.tsv",
            "AMR_DNA-Klebsiella_oxytoca.tsv",
            "AMR_DNA-Klebsiella_pneumoniae.tsv",
            "AMR_DNA-Neisseria_gonorrhoeae.tsv",
            "AMR_DNA-Salmonella.tsv",
            "AMR_DNA-Staphylococcus_aureus.tsv",
            "AMR_DNA-Streptococcus_pneumoniae.tsv",
        ]
        dna_fas = [
            "AMR_DNA-Acinetobacter_baumannii.fa",
            "AMR_DNA-Campylobacter.fa",
            "AMR_DNA-Clostridioides_difficile.fa",
            "AMR_DNA-Enterococcus_faecalis.fa",
            "AMR_DNA-Enterococcus_faecium.fa",
            "AMR_DNA-Escherichia.fa",
            "AMR_DNA-Klebsiella_oxytoca.fa",
            "AMR_DNA-Klebsiella_pneumoniae.fa",
            "AMR_DNA-Neisseria_gonorrhoeae.fa",
            "AMR_DNA-Salmonella.fa",
            "AMR_DNA-Staphylococcus_aureus.fa",
            "AMR_DNA-Streptococcus_pneumoniae.fa",
        ]
        for tsv in dna_tsvs:
            sh(f"wget -nv -c {base}/{tsv} 1>> databaseprep.log 2>> databaseprep.error")
        for fa in dna_fas:
            sh(f"wget -nv -c {base}/{fa} 1>> databaseprep.log 2>> databaseprep.error")
        # cat and drop header lines with 'mutation_position'
        with open("AMR_DNA-mutation.tsv", "w", encoding="utf-8") as out:
            for tsv in dna_tsvs:
                with open(tsv, "r", encoding="utf-8") as fin:
                    for ln in fin:
                        if "mutation_position" in ln:
                            continue
                        out.write(ln)
        # cat fas lines
        with open("AMR_DNA.fa", "w", encoding="utf-8") as out:
            for fa in dna_fas:
                with open(fa, "r", encoding="utf-8") as fin:
                    for ln in fin:
                        out.write(ln)

    # Build AMR_DNA_underscore_combined.tsv
    if not exists("AMR_DNA_underscore_combined.tsv"):
        print("Prepping AMR_DNA mutation database")
        build_AMR_DNA_combined(database)

    # Build KMA DB
    if not exists("AMR_DNA_underscore.comp.b"):
        print("Building AMR_DNA KMA database")
        sh("kma index -i AMR_DNA_underscore.fa -o AMR_DNA_underscore 1>> databaseprep.log 2>> databaseprep.error")

    # Populate the global class set so summaries include empty classes
    global ALL_CLASSES
    ALL_CLASSES = collect_all_classes(database)



    # ------------------------
    # Alignments
    # ------------------------
    os.chdir(output)
    sample_fastq = output / f"{sample}.fastq"

    # Prepare some initial files (class/accession lists)
    print("Preparing initial class/accession files")
    
    # AMRProt classes (col 6) excluding header "class"
    with open(f"{database}/AMRProt-mutation_underscore.tsv", encoding="utf-8") as f:
        classes = {line.rstrip("\n").split("\t")[5] for line in f if not line.startswith("class")}
    with open("class", "w", encoding="utf-8") as fout:
        for c in sorted(classes):
            fout.write(c + "\n")
    
    # AMRProt accessions (col 2)
    with open(f"{database}/AMRProt-mutation_underscore.tsv", encoding="utf-8") as f:
        accs = {line.rstrip("\n").split("\t")[1] for line in f if line.strip()}
    with open("accession", "w", encoding="utf-8") as fout:
        for a in sorted(accs):
            fout.write(a + "\n")
    
    # AMR_DNA classes (col 5) excluding header "class"
    with open(f"{database}/AMR_DNA-mutation_underscore.tsv", encoding="utf-8") as f:
        dclasses = {line.rstrip("\n").split("\t")[4] for line in f if not line.startswith("class")}
    with open("dna_class", "w", encoding="utf-8") as fout:
        for dc in sorted(dclasses):
            fout.write(dc + "\n")
    
    # AMR_DNA accessions (col 1)
    with open(f"{database}/AMR_DNA-mutation_underscore.tsv", encoding="utf-8") as f:
        daccs = {line.rstrip("\n").split("\t")[0] for line in f if line.strip()}
    with open("dna_accession", "w", encoding="utf-8") as fout:
        for da in sorted(daccs):
            fout.write(da + "\n")
            
    # DIAMOND
    print("Aligning reads to proteins using diamond blastx")
    sh(
        f"diamond blastx --threads {threads} -d {database}/AMRProt -q {sample_fastq} "
        f"-o {sample}.prot.hits.txt -F 15 --range-culling -k1 --range-cover 5 "
        f"--id {identity} --masking 0 --outfmt 6 "
        f"qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore full_sseq qseq_translated "
        f"1>> {sample}.log 2>> {sample}.error"
    )

    # Build protein input (exact accession match via dict)
    print("Parsing diamond output")

    out_path = f"{sample}.prot.input.tsv"
    with open(out_path, "w", encoding="utf-8") as fout:
        fout.write("class\tgene\tread\treference\ttarget\tchanges_str\n")

    # Index AMRProt combined: accession -> list[(class, changes)]
    amr_map = {}
    with open(f"{database}/AMRProt-mutation_underscore_combined.tsv", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                cls, acc, chg = parts[0], parts[1], parts[2]
                amr_map.setdefault(acc, []).append((cls, chg))

    def pick_fields_from_hit(line: str):
        # Get the right columns
        cols = line.rstrip("\n").split("\t")
        tmp = [(cols[i] if i < len(cols) else "") for i in (0, 1, 12, 13)]
        expanded = "\t".join(tmp).replace("|", "\t").split("\t")
        if len(expanded) < 14:
            expanded += [""] * (14 - len(expanded))
        _id        = expanded[0]
        _accession = expanded[2]
        _gene      = expanded[11]
        _reference = expanded[12]
        _target    = expanded[13]
        return _id, _accession, _gene, _reference, _target

    with open(f"{sample}.prot.hits.txt", encoding="utf-8") as fh, open(out_path, "a", encoding="utf-8") as fout:
        for line in fh:
            _id, _acc, _gene, _ref, _tgt = pick_fields_from_hit(line)
            if not _acc:
                continue
            for cls, chg in amr_map.get(_acc, []):
                fout.write("\t".join([cls, _gene, _id, _ref, _tgt, chg]) + "\n")



    # Score proteins (R)
    print("Scoring amino acid substitutions")
    #shutil.copyfile(f"{sample}.prot.input.tsv", "input.tsv")
    #sh(f'R --vanilla < "{script_dir}/prot_score_mutations.R" 1>> {sample}.log 2>> {sample}.error')
    sh(f'Rscript "{script_dir}/prot_score_mutations.R" "{sample}.prot.input.tsv" {threads} 1>> {sample}.log 2>> {sample}.error')
    shutil.copyfile("updated_table_with_scores_and_mutations.tsv", f"{sample}.prot.updated_table_with_scores_and_mutations.tsv")
    #os.remove("input.tsv")
    os.remove("updated_table_with_scores_and_mutations.tsv")

    # KMA
    print("Aligning reads to genes using kma")
    identityfloat = float(identity) / 100
    sh(
        f'kma -bcNano -hmm -ont -t_db {database}/AMR_DNA_underscore -i {sample_fastq} '
        f'-o {sample} -t {threads} -nc -na -1t1 -mrs {identityfloat} 1>> {sample}.log 2>> {sample}.error'
    )
    # delete temp fastq
    os.remove(sample_fastq)

    # DNA input 
    
    print("Parsing KMA output")
    out_path = f"{sample}.dna.input.tsv"
    with open(out_path, "w", encoding="utf-8") as fout:
        fout.write("class\tgene\tread\treference\ttarget\tchanges_str\n")
    
    # Index AMR_DNA combined table: accession -> list[(class, reference, changes)]
    dna_map = {}
    with open(f"{database}/AMR_DNA_underscore_combined.tsv", encoding="utf-8") as fdb:
        for line in fdb:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                cls, acc, ref, chg = parts[0], parts[1], parts[2], parts[3]
                dna_map.setdefault(acc, []).append((cls, ref, chg))
    
    # Read KMA fragments (.frag.gz)
    # take columns 1,6,7 (0-based: 0,5,6); trim the read at first space.
    with gzip.open(f"{sample}.frag.gz", "rt", encoding="utf-8") as ffrag, open(out_path, "a", encoding="utf-8") as fout:
        for line in ffrag:
            cols = line.rstrip("\n").split("\t")
            # Ensure we have at least 7 columns
            if len(cols) < 7:
                continue
            target = cols[0]
            accession = cols[5]
            read_id = cols[6].split(" ")[0]  # trim at first space
    
            if not accession:
                continue
    
            for cls, ref, chg in dna_map.get(accession, []):
                # NOTE: "gene" column for DNA is the accession
                fout.write("\t".join([cls, accession, read_id, ref, target, chg]) + "\n")
    



    # Score DNA (R)
    print("Scoring DNA mutations")
    #shutil.copyfile(f"{sample}.dna.input.tsv", "input.tsv")
    #sh(f'R --vanilla < "{script_dir}/dna_score_mutations.R" 1>> {sample}.log 2>> {sample}.error')
    sh(f'Rscript "{script_dir}/dna_score_mutations.R" "{sample}.dna.input.tsv" {threads} 1>> {sample}.log 2>> {sample}.error')
    shutil.copyfile("updated_table_with_scores_and_mutations.tsv", f"{sample}.dna.updated_table_with_scores_and_mutations.tsv")
    os.remove("updated_table_with_scores_and_mutations.tsv")

    # ------------------------
    # Summaries
    # ------------------------
    write_class_and_gene_summaries(sample, output, is_prot=True)
    write_class_and_gene_summaries(sample, output, is_prot=False)

    # Reporting
    print()
    if exists(f"{sample}.prot.updated_table_with_scores_and_mutations.tsv"):
        print(f"Protein read classification output (WT/R) can be found in the file {output}/{sample}.prot.updated_table_with_scores_and_mutations.tsv")
    else:
        print("Protein read classification output not generated. Something went wrong")

    if exists(f"{sample}.class.prot.summary.txt"):
        print(f"Class protein summary of classification output (WT/R) can be found in the file {output}/{sample}.class.prot.summary.txt")
        print(f"Gene protein summary of classification output (WT/R) can be found in the file {output}/{sample}.gene.prot.summary.txt")
        print()
    else:
        print("Protein summary output not generated. Something went wrong")

    if exists(f"{sample}.dna.updated_table_with_scores_and_mutations.tsv"):
        print(f"DNA read classification output (WT/R) can be found in the file {output}/{sample}.dna.updated_table_with_scores_and_mutations.tsv")
    else:
        print("DNA read classification output not generated. Something went wrong")

    if exists(f"{sample}.class.dna.summary.txt"):
        print(f"Class DNA summary of classification output (WT/R) can be found in the file {output}/{sample}.class.dna.summary.txt")
        print(f"Gene DNA summary of classification output (WT/R) can be found in the file {output}/{sample}.gene.dna.summary.txt")
        print()
        sys.exit(0)
    else:
        print("DNA summary output not generated. Something went wrong")

if __name__ == "__main__":
    main()
