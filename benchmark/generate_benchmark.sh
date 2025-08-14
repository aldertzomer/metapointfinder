# this generates a set of benchmark fastq files based on the amrprot database of amrfinder as processed by metapointfinder
# copy the files in this folder to the metapointfinderdb folder and run it from there

#make output folder
mkdir benchmark

#generate 3 randomly selected mutated protein fasta files from the amrprot database. also include 3 wt sequences
python  mutation_r_wt_generator.py AMRProt-mutation_underscore_combined.tsv AMRProt_mutation.fa benchmark/mutated.faa --per-row 3 --seed 42

# make it a bit simpler
cat benchmark/mutated.faa |sed 's/|/_/g' |sed 's/=/_/g' > benchmark/mutated_underscore.faa

# backtranslate the protein fasta to dna using a variety of codon usage tables
for cut in Eacica Eagrtu Ebacsu Ebraja Ecaucr Echltr Echlre Ecloab Eerwct Eecoli Ebacst Ehaein Ehalsa Eklepn Elacdl Emyctu Eneigo Epseae Erhile Erhosh Esalty Eserma Erhime Estaau Estrpn Estrco Esynco Etheth Evibch Eyeren  ; do
    echo "backtranseq -seq benchmark/mutated_underscore.faa -outfile benchmark/$cut.fasta -cfile /home/zomer007/miniforge3/envs/genomics/share/EMBOSS/data/CODONS/$cut.cut"
done |parallel -j 64

# pad the sequences to 5 kb to simulate 5 kb pseudogenomes
for cut in Eacica Eagrtu Ebacsu Ebraja Ecaucr Echltr Echlre Ecloab Eerwct Eecoli Ebacst Ehaein Ehalsa Eklepn Elacdl Emyctu Eneigo Epseae Erhile Erhosh Esalty Eserma Erhime Estaau Estrpn Estrco Esynco Etheth Evibch Eyeren  ; do
    echo "python pad_to_10kb.py benchmark/$cut.fasta benchmark/$cut.5001.fasta --length 5001 --trim center --seed 42"
done |parallel -j 64

# generate reads with different lengths using wgsim with 1% error rate from the 5 kb pseudogenomes
for cut in Eacica Eagrtu Ebacsu Ebraja Ecaucr Echltr Echlre Ecloab Eerwct Eecoli Ebacst Ehaein Ehalsa Eklepn Elacdl Emyctu Eneigo Epseae Erhile Erhosh Esalty Eserma Erhime Estaau Estrpn Estrco Esynco Etheth Evibch Eyeren  ; do
    echo "wgsim -N 100000 -1 100 -d 0 -S 42 -e 0 -r 0.01 benchmark/$cut.5001.fasta benchmark/$cut.100.fastq /dev/null"
    echo "wgsim -N 50000 -1 200 -d 0 -S 42 -e 0 -r 0.01 benchmark/$cut.5001.fasta benchmark/$cut.200.fastq /dev/null"
    echo "wgsim -N 10000 -1 1000 -d 0 -S 42 -e 0 -r 0.01 benchmark/$cut.5001.fasta benchmark/$cut.1000.fastq /dev/null"
    echo "wgsim -N 2000 -1 5000 -d 0 -S 42 -e 0 -r 0.01 benchmark/$cut.5001.fasta benchmark/$cut.5000.fastq /dev/null"
done |parallel -j 64 >/dev/null


