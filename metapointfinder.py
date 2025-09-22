#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
metapointfinder.py
Reimplementation of the Bash metapointfinder pipeline in Python, preserving
all quirks of formatting, sorting, and mutation handling.

Usage:
  python metapointfinder.py <file.fastq[.gz]> <databasefolder> <outputfolder> [force]

Requirements:
  - Programs listed in ./dependencies (diamond, kma, wget, R, zcat, cut, sed, grep, tr, sort, uniq, awk, tail, gzip)
  - R packages: msa, Biostrings, pwalign
  - DIAMOND 2.0.15 (as in your README requirement)

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
    Accession is the header up to first whitespace (matches Bash behavior).
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
# Mutation combiners (Python versions of the Bash loops)
# ------------------------

def parse_amrprot_mutation_underscore_combined(tsv_in):
    """
    Build dict: (class, accession) -> sorted, comma-joined changes_str (low->high).
    Behavior mirrors the Bash:
      - Use the 'position' column (numeric) from AMRProt-mutation.tsv (cut -f 2,3,4,6)
      - Keep raw token until the END, then map 'del'/'STOP' -> '-'
      - Do NOT strip '-' from the right side (so '...107-' is preserved)
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
    Mimic the Bash AMR_DNA combiner exactly:
      - Replace '/' by '_' in mutation TSV and FASTA headers
      - Build AMR_DNA_underscore_fa.tsv (accession \t sequence)
      - For each (class, accession), accumulate mutations, build changes_str
        (low->high) and only THEN map 'del'/'STOP' to '-'
      - Do NOT remove '-' from right side so trailing '-' deletions are preserved
      - Sort output lines by (class, accession)
    """
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
    rows = []
    with open(dna_mut_us, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            if "mutation_position" in line or ("class" in line and "position" in line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            acc  = parts[0]
            pos  = parts[1]
            mut  = parts[2]
            cls  = parts[4]
            rows.append((acc, pos, mut, cls))

    classes    = sorted(set(cls for *_, cls in rows if cls != "class"))
    accessions = sorted(set(acc for acc, *_ in rows))

    by_cls_acc = {}
    for cls in classes:
        for acc in accessions:
            muts = []
            for (acc_i, pos, mut, cls_i) in rows:
                if acc_i != acc or cls_i != cls:
                    continue

                mutation_cleaned = mut.split("_", 1)[1] if "_" in mut else mut

                mpos = re.search(r"[0-9]+", pos)
                if not mpos:
                    continue
                pos_int = int(mpos.group(0))

                pieces = re.split(r"([0-9]+)", mutation_cleaned, maxsplit=1)
                if len(pieces) < 3:
                    continue
                baseleft  = pieces[0]
                baseright = pieces[2]  # keep any '-'

                mutation_raw = f"{baseleft}{pos_int}{baseright}"
                muts.append((pos_int, mutation_raw))

            if muts:
                muts.sort(key=lambda x: x[0])
                changes_str = ",".join(m for _, m in muts)
                # Only now map del/STOP
                changes_str = changes_str.replace("del", "-").replace("STOP", "-")
                seq = acc_to_seq.get(acc, "")
                by_cls_acc[(cls, acc)] = (seq, changes_str)

    # Write sorted by (class, accession)
    with open(combined_out, "w", encoding="utf-8") as out:
        for (cls, acc) in sorted(by_cls_acc.keys(), key=lambda t: (t[0], t[1])):
            seq, changes = by_cls_acc[(cls, acc)]
            out.write(f"{cls}\t{acc}\t{seq}\t{changes}\n")

# ------------------------
# Summaries (Python version of the Bash counts)
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
    classes = sorted(set(r["class"] for r in rows if r.get("class") and r["class"] != "class"))
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
    print("This is metapointfinder v0.2")
    print("This tool finds substitions in translated reads matching to relevant proteins and mutations in reads matching to relevant genes")
    print("Relevant proteins and genes are obtained from the AMRFinder database")

    if len(sys.argv) < 4:
        print("usage: metapointfinder.py file.fastq[.gz] databasefolder outputfolder [force]")
        sys.exit(1)

    script_path = Path(__file__).resolve()
    script_dir  = script_path.parent

    input_path  = Path(sys.argv[1]).resolve()
    database    = Path(sys.argv[2])
    output      = Path(sys.argv[3])
    force_flag  = sys.argv[4] if len(sys.argv) > 4 else None

    fastqfile = input_path.name

    # ------------------------
    # Dependencies (from ./dependencies)
    # ------------------------
    deps_file = script_dir / "dependencies"
    if not deps_file.exists():
        print("Missing 'dependencies' file next to the script.")
        sys.exit(1)

    # Check each dependency (like Bash)
    with open(deps_file, "r", encoding="utf-8") as f:
        for program in [ln.strip() for ln in f if ln.strip()]:
            if shutil.which(program) is None:
                print(f"Required dependency {program} not installed")
                sys.exit(1)

    # R packages
    try:
        sh('R -e "stopifnot(3 == length(find.package(c(\'msa\', \'Biostrings\', \'pwalign\'))))" --slave')
    except subprocess.CalledProcessError:
        print("R package msa and/or Biostrings and/or pwalign not installed")
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

    # AMRProt
    if not exists("AMRProt.fa"):
        print("downloading https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa")
        sh("wget -nv -c https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa 1>> databaseprep.log 2>> databaseprep.error")

    if not exists("AMRProt-mutation.tsv"):
        print("downloading https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt-mutation.tsv")
        sh("wget -nv -c https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt-mutation.tsv 1>> databaseprep.log 2>> databaseprep.error")

    # Prepare AMRProt-mutation_underscore.tsv (slashes -> underscores)
    if not exists("AMRProt-mutation_underscore.tsv"):
        s = read_text("AMRProt-mutation.tsv").replace("/", "_")
        write_text("AMRProt-mutation_underscore.tsv", s)

    # Combined AMRProt mutations (Python replicant of the Bash loop)
    if not exists("AMRProt-mutation_underscore_combined.tsv"):
        print("Prepping AMRProt mutation database")
        combined_map = parse_amrprot_mutation_underscore_combined("AMRProt-mutation.tsv")
        # write sorted by (class, accession)
        with open("AMRProt-mutation_underscore_combined.tsv", "w", encoding="utf-8") as out:
            for (cls, acc) in sorted(combined_map.keys(), key=lambda t: (t[0], t[1])):
                out.write(f"{cls}\t{acc}\t{combined_map[(cls, acc)]}\n")

    # Build AMRProt_mutation.fa (only headers containing "mutation")
    if not exists("AMRProt_mutation.fa"):
        # Equivalent of: cat AMRProt.fa |tr "\n" "\t" |tr ">" "\n" |grep mutation |sed 's/^/>/'|tr "\t" "\n"
        lines = []
        with open("AMRProt.fa", "r", encoding="utf-8") as f:
            buf = f.read()
        # Split records
        recs = [r for r in buf.split(">") if r]
        for rec in recs:
            if "mutation" in rec:
                lines.append(">" + rec)
        write_text("AMRProt_mutation.fa", "".join(lines))

    # Build diamond DB
    if not exists("AMRProt.dmnd"):
        print("Building AMRProt diamond database")
        sh("diamond makedb --in AMRProt_mutation.fa --db AMRProt 1>> databaseprep.log 2>> databaseprep.error")

    # DNA databases
    if not exists("AMR_DNA.fa"):
        print("downloading DNA fasta databases")
        base = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
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
        for fna in dna_fas:
            sh(f"wget -nv -c {base}/{fna} 1>> databaseprep.log 2>> databaseprep.error")
        # concat
        with open("AMR_DNA.fa", "wb") as fout:
            for fna in dna_fas:
                with open(fna, "rb") as fin:
                    shutil.copyfileobj(fin, fout)

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
        for tsv in dna_tsvs:
            sh(f"wget -nv -c {base}/{tsv} 1>> databaseprep.log 2>> databaseprep.error")
        # cat and drop header lines with 'mutation_position'
        with open("AMR_DNA-mutation.tsv", "w", encoding="utf-8") as out:
            for tsv in dna_tsvs:
                with open(tsv, "r", encoding="utf-8") as fin:
                    for ln in fin:
                        if "mutation_position" in ln:
                            continue
                        out.write(ln)

    # Build AMR_DNA_underscore_combined.tsv (Python exact replica of Bash loop)
    if not exists("AMR_DNA_underscore_combined.tsv"):
        print("Prepping AMR_DNA mutation database")
        build_AMR_DNA_combined(database)

    # Build KMA DB
    if not exists("AMR_DNA_underscore.comp.b"):
        print("Building AMR_DNA KMA database")
        sh("kma index -i AMR_DNA_underscore.fa -o AMR_DNA_underscore 1>> databaseprep.log 2>> databaseprep.error")

    # ------------------------
    # Alignments
    # ------------------------
    os.chdir(output)
    sample_fastq = output / f"{sample}.fastq"

    # Prepare some initial files (class/accession lists)
    # This mirrors the Bash cuts; we’ll reuse only where needed for summaries later
    sh(f"cat {database}/AMRProt-mutation_underscore.tsv |cut -f 6 |sort |uniq |grep -v class > class")
    sh(f"cat {database}/AMRProt-mutation_underscore.tsv |cut -f 2 |sort |uniq > accession")
    sh(f"cat {database}/AMR_DNA-mutation_underscore.tsv |cut -f 5 |sort |uniq |grep -v class > dna_class")
    sh(f"cat {database}/AMR_DNA-mutation_underscore.tsv |cut -f 1 |sort |uniq > dna_accession")

    # DIAMOND
    print("Aligning reads to proteins using diamond blastx")
    sh(
        f"diamond blastx -d {database}/AMRProt -q {sample_fastq} "
        f"-o {sample}.prot.hits.txt -F 15 --range-culling -k1 --range-cover 5 --iterate "
        f"--id 70 --masking 0 --outfmt 6 "
        f"qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore full_sseq qseq_translated "
        f"1>> {sample}.log 2>> {sample}.error"
    )

    # Build protein input (REPLICATE the exact shell pipeline to preserve field handling)
    print("Parsing diamond output")
    with open(f"{sample}.prot.input.tsv", "w", encoding="utf-8") as fout:
        fout.write("class\tgene\tread\treference\ttarget\tchanges_str\n")
    sh(
        f"""cat {sample}.prot.hits.txt |cut -f 1,2,13,14 |sed "s/|/\\t/g" |cut -f 1,3,12,13,14 | \
while read id accession gene reference target ; do \
  cat {database}/AMRProt-mutation_underscore_combined.tsv |grep $accession |while read class accession changes ; do \
    echo $class $gene $id $reference $target $changes ; \
  done ; \
done | tr " " "\\t" >> {sample}.prot.input.tsv"""
    )

    # Score proteins (R)
    print("Scoring amino acid substitutions")
    shutil.copyfile(f"{sample}.prot.input.tsv", "input.tsv")
    sh(f'R --vanilla < "{script_dir}/prot_score_mutations.R" 1>> {sample}.log 2>> {sample}.error')
    shutil.copyfile("updated_table_with_scores_and_mutations.tsv", f"{sample}.prot.updated_table_with_scores_and_mutations.tsv")
    os.remove("input.tsv")
    os.remove("updated_table_with_scores_and_mutations.tsv")

    # KMA
    print("Aligning reads to genes using kma")
    sh(
        f'kma -t -bcNano -hmm -ont -t_db {database}/AMR_DNA_underscore -i {sample_fastq} '
        f'-o {sample} -t 16 -nc -na -1t1 1>> {sample}.log 2>> {sample}.error'
    )
    # delete temp fastq
    os.remove(sample_fastq)

    # DNA input (replicate Bash pipeline on frag.gz)
    print("Parsing KMA output")
    with open(f"{sample}.dna.input.tsv", "w", encoding="utf-8") as fout:
        fout.write("class\tgene\tread\treference\ttarget\tchanges_str\n")
    sh(
        f"""zcat {sample}.frag.gz |cut -f 1,6,7 |cut -f 1 -d " " | \
while read target accession read ; do \
  cat {database}/AMR_DNA_underscore_combined.tsv |grep $accession |while read class accession reference changes; do \
    echo $class $accession $read $reference $target $changes ; \
  done ; \
done | tr " " "\\t" >> {sample}.dna.input.tsv"""
    )

    # Score DNA (R)
    print("Scoring DNA mutations")
    shutil.copyfile(f"{sample}.dna.input.tsv", "input.tsv")
    sh(f'R --vanilla < "{script_dir}/dna_score_mutations.R" 1>> {sample}.log 2>> {sample}.error')
    shutil.copyfile("updated_table_with_scores_and_mutations.tsv", f"{sample}.dna.updated_table_with_scores_and_mutations.tsv")
    os.remove("input.tsv")
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
