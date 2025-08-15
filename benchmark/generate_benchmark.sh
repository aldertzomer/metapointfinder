# this generates a set of benchmark fastq files based on the amrprot database of amrfinder as processed by metapointfinder
# copy the files in this folder to the metapointfinderdb folder and run it from there

#make output folder
mkdir benchmark

#generate 3 randomly selected mutated protein fasta files from the amrprot database. also include 3 wt sequences
python  mutation_r_wt_generator.py AMRProt-mutation_underscore_combined.tsv AMRProt_mutation.fa benchmark/mutated.faa --per-row 3 --seed 42

# make it a bit simpler
cat benchmark/mutated.faa |sed 's/|/_/g' |sed 's/=/_/g' > benchmark/mutated_underscore.faa

# backtranslate the protein fasta to dna using a variety of codon usage tables
#for cut in Eacica Eagrtu Ebacsu Ebraja Ecaucr Echltr Echlre Ecloab Eerwct Eecoli Ebacst Ehaein Ehalsa Eklepn Elacdl Emyctu Eneigo Epseae Erhile Erhosh Esalty Eserma Erhime Estaau Estrpn Estrco Esynco Etheth Evibch Eyeren  ; do
for cut in Eecoli Estaau Emyctu ; do
    echo "backtranseq -seq benchmark/mutated_underscore.faa -outfile benchmark/$cut.fasta -cfile /home/zomer007/miniforge3/envs/genomics/share/EMBOSS/data/CODONS/$cut.cut"
done |parallel -j 64

# pad the sequences to 6 kb to simulate 6 kb pseudogenomes
#for cut in Eacica Eagrtu Ebacsu Ebraja Ecaucr Echltr Echlre Ecloab Eerwct Eecoli Ebacst Ehaein Ehalsa Eklepn Elacdl Emyctu Eneigo Epseae Erhile Erhosh Esalty Eserma Erhime Estaau Estrpn Estrco Esynco Etheth Evibch Eyeren  ; do
for cut in Eecoli Estaau Emyctu ; do
    echo "python pad_to_10kb.py benchmark/$cut.fasta benchmark/$cut.6000.fasta --length 6000 --trim center --seed 42"
done |parallel -j 64

# generate reads with different lengths using wgsim with varying error rates from the 6 kb pseudogenomes
#for cut in Eacica Eagrtu Ebacsu Ebraja Ecaucr Echltr Echlre Ecloab Eerwct Eecoli Ebacst Ehaein Ehalsa Eklepn Elacdl Emyctu Eneigo Epseae Erhile Erhosh Esalty Eserma Erhime Estaau Estrpn Estrco Esynco Etheth Evibch Eyeren  ; do
for cut in Eecoli Estaau Emyctu ; do
    for r in 0.00 0.01 0.02 0.03 0.05 0.1 0.2 ; do
	echo "wgsim -N 100000 -1 100 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.100.$r.fastq /dev/null"
	echo "wgsim -N 50000 -1 200 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.200.$r.fastq /dev/null"
	echo "wgsim -N 10000 -1 1000 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.1000.$r.fastq /dev/null"
	echo "wgsim -N 2000 -1 5000 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.5000.$r.fastq /dev/null"
    done
done |parallel -j 64 >/dev/null

# combine the data from the different codon usage tables into one file
for r in 0.00 0.01 0.02 0.03 0.05 0.1 0.2 ; do
    cat benchmark/E*.100.$r.fastq > benchmark/allspecies.100.$r.fastq
    cat benchmark/E*.200.$r.fastq > benchmark/allspecies.200.$r.fastq
    cat benchmark/E*.1000.$r.fastq > benchmark/allspecies.1000.$r.fastq
    cat benchmark/E*.5000.$r.fastq > benchmark/allspecies.5000.$r.fastq
done
