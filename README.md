# MetaPointFinder

**MetaPointFinder** is a tool for detecting and scoring resistance-associated point mutations directly from long-read sequencing data (e.g. nanopore reads), using the AMRFinder database as reference. It combines protein- and DNA-level mutation detection to classify reads as wild-type (WT) or resistant (R), providing both read-level and summary-level output.

---

## Features

- Detects amino acid substitutions in translated reads using DIAMOND and AMRFinder protein mutation databases.
- Detects nucleotide mutations using KMA against AMRFinder DNA databases.
- Scores mutations using MSA-based alignment and outputs both read-level classifications and summary tables.
- Supports direct analysis of long-read metagenomic datasets (nanopore fastq files).

---

## Installation

### Dependencies

- **R** with libraries:
  - `msa`
  - `Biostrings`

- **DIAMOND v2.0.15** (this exact version is required)
- **KMA**
- **wget**
- Standard Unix utilities (`awk`, `sed`, `cut`, etc.)

R packages can be installed with:

```r
install.packages("BiocManager")
BiocManager::install(c("msa", "Biostrings"))
```

DIAMOND v2.0.15 can be downloaded from:
https://github.com/bbuchfink/diamond/releases/tag/v2.0.15

KMA can be obtained from:
https://bitbucket.org/genomicepidemiology/kma

### Installation

```bash
git clone https://github.com/aldertzomer/metapointfinder.git
cd metapointfinder
chmod +x metapointfinder
```

---

## Usage

```bash
./metapointfinder input.fastq database_folder output_folder
```

OR for gzipped input:

```bash
./metapointfinder input.fastq.gz database_folder output_folder
```

### Arguments

- `input.fastq` or `input.fastq.gz`: Nanopore long-read data file.
- `database_folder`: Directory where AMRFinder databases will be downloaded/prepared.
- `output_folder`: Directory to store all results.

> If `database_folder` does not exist, it will be created automatically and the required AMRFinder databases will be downloaded and processed.

---

## Pipeline Overview

1. **Database Preparation:**
   - AMRFinder protein and DNA reference databases are downloaded.
   - Mutation tables are preprocessed for scoring.
   - DIAMOND and KMA databases are created.

2. **Read Processing:**
   - Reads aligned to AMR proteins using DIAMOND blastx.
   - Reads aligned to AMR genes using KMA.
   - Matching regions extracted and scored using the provided R scripts.

3. **Mutation Scoring (R):**
   - Translated protein reads and DNA reads are aligned to references using ClustalW (via the `msa` R package).
   - Known mutations from AMRFinder mutation tables are searched for.
   - Each read is scored:
     - **MutationScore:** Number of known resistance mutations detected.
     - **DetectedMutations:** List of mutations found (or "None").

4. **Results Aggregation:**
   - WT and R counts summarized by antibiotic class and gene.

---

## Output Files

| File                                               | Description                                                     |
|----------------------------------------------------|-----------------------------------------------------------------|
| `*.prot.updated_table_with_scores_and_mutations.tsv` | Read-level protein mutation classification (WT/R).              |
| `*.dna.updated_table_with_scores_and_mutations.tsv`  | Read-level DNA mutation classification (WT/R).                  |
| `*.class.prot.summary.txt`                           | WT/R summary per antibiotic class (protein level).              |
| `*.gene.prot.summary.txt`                            | WT/R summary per gene (protein level).                          |
| `*.class.dna.summary.txt`                            | WT/R summary per antibiotic class (DNA level).                  |
| `*.gene.dna.summary.txt`                             | WT/R summary per gene (DNA level).                              |
| Log and error files                                  | Detailed logs for troubleshooting DIAMOND, KMA, and R execution.|

---

## Example Output (Class Summary)

Example (`*.class.prot.summary.txt`):

```
class           WT    R
BETA-LACTAM     255   44
COLISTIN        59    12
QUINOLONE       185   25
MULTIDRUG       134   33
```

---

## Example Output (Read Classification)

Example (`*.prot.updated_table_with_scores_and_mutations.tsv`):

```
class\tgene\tread\treference\ttarget\tchanges_str\tMutationScore\tDetectedMutations
BETA-LACTAM\ttwo-component_system_sensor_histidine_kinase_BaeS\t10478f68-f181-4626-976f-93d14c49844b\t...\t...\tY42H,T175P,R153P\t0\tNone
COLISTIN\ttwo-component_system_sensor_histidine_kinase_PmrB\t1249bd62-efe5-46d5-96c4-1e903c85dec5\t...\t...\tV161G,T92P\t1\tV161G
...
```

- **class:** Antibiotic class.
- **gene:** Gene or protein name.
- **read:** Read identifier.
- **reference:** Reference sequence used for alignment.
- **target:** Sequence from read.
- **changes_str:** List of expected mutations.
- **MutationScore:** Number of detected resistance mutations.
- **DetectedMutations:** Mutations actually found in the read.

---

## Notes

- Ensure **DIAMOND v2.0.15** is used. Other versions may not work as expected.
- Databases are downloaded directly from NCBI during the first run.
- Mutation detection relies on known mutations in the AMRFinder database. Novel mutations will not be scored.
- For metagenomic samples with highly fragmented reads, alignment quality directly affects detection sensitivity.

---

## License

This tool is released under the Apache-2.0 License.

---

## Citation

If you use MetaPointFinder in your work, please cite:

> Zomer, A.L. *et al.* MetaPointFinder: a tool to detect resistance mutations directly from long-read metagenomic data using AMRFinder databases. [GitHub Repository](https://github.com/aldertzomer/metapointfinder)
