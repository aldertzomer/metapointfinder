#!/bin/bash
# this generates a set of benchmark fastq files based on the amr_dna database of amrfinder as processed by metapointfinder and then processes these and scores them for TP, TN, FP, etc
# copy the files in this folder to the metapointfinderdb folder and run it from there

# set location of tool and db
tool=/mnt/data/tools/metapointfinder/metapointfinder/metapointfinder.py
db=/mnt/data/db/metapointfinderdb

# make output dir
mkdir benchmarkdna

# generate mutant and wt sequences
python make_dna_mutants.py AMR_DNA_underscore_combined.tsv benchmarkdna/mutated.fna --per-row 3 --seed 42

# make it a bit simpler
cat benchmarkdna/mutated.fna |sed 's/|/_/g' |sed 's/=/_/g' > benchmarkdna/mutated_underscore.fna

# pad the sequences to 6 kb to simulate 6 kb pseudogenomes
python pad_to_10kb.py benchmarkdna/mutated_underscore.fna benchmarkdna/mutated_underscore.6000.fasta --length 6000 --trim center --seed 42

# generate error prone reads of different lengths
for r in 0.00 0.01 0.02 0.03 0.05 0.1 0.2 0.3 ; do
    echo "wgsim -N 1000000 -1 100 -d 0 -S 42 -e 0 -r $r benchmarkdna/mutated_underscore.6000.fasta benchmarkdna/100.$r.fastq /dev/null"
    echo "wgsim -N 500000 -1 200 -d 0 -S 42 -e 0 -r $r benchmarkdna/mutated_underscore.6000.fasta benchmarkdna/200.$r.fastq /dev/null"
    echo "wgsim -N 100000 -1 1000 -d 0 -S 42 -e 0 -r $r benchmarkdna/mutated_underscore.6000.fasta benchmarkdna/1000.$r.fastq /dev/null"
    echo "wgsim -N 20000 -1 5000 -d 0 -S 42 -e 0 -r $r benchmarkdna/mutated_underscore.6000.fasta benchmarkdna/5000.$r.fastq /dev/null"
done |parallel -j 64 >/dev/null

# run metapointfinder on all outputs. change the locations of the tool and the db to suit your needs
for r in 0.00 0.01 0.02 0.03 0.05 0.1 0.2 0.3 ; do
        echo "python $tool --input benchmarkdna/100.$r.fastq --db $db --output benchmarkdna/mpf.100.$r --identity 50 --threads 2"
        echo "python $tool --input benchmarkdna/200.$r.fastq --db $db --output benchmarkdna/mpf.200.$r --identity 50 --threads 2"
        echo "python $tool --input benchmarkdna/1000.$r.fastq --db $db --output benchmarkdna/mpf.1000.$r --identity 50 --threads 2"
	echo "python $tool --input benchmarkdna/5000.$r.fastq --db $db --output benchmarkdna/mpf.5000.$r --identity 50 --threads 2"
done |parallel -j 64

# this thing checks if the information in the header is the same as the prediction. this counts unique hits. We only check WT reads and Wildtype prediction and RES reads and Resistance prediction. There will be many false negatives and unknowns as many of the shorter reads do not cover the resistance determining region, therefore these are not scored. False Positives can occur by chance in the reads with mutation.
echo "file,WT,RES,FP,FN,WT_UNK,RES_UNK" > benchmarkdna/final_scores.csv
find benchmarkdna |grep dna.updated_table_with_scores_and_mutations.tsv$ |while read output ; do
    echo $output |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Wildtype |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Wildtype |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Unknown |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Unknown |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l
done >>benchmarkdna/final_scores.csv

# this thing checks if the information in the header is the same as the prediction. This counts all hits. We only check WT reads and Wildtype prediction (TN) and RES reads and Resistance prediction (TP). There will be many false negatives and unknowns as many of the shorter reads do not cover the resistance determining region, therefore these are not scored. False Positives can occur by chance in the reads with mutation.
echo "file,WT,RES,FP,FN,WT_UNK,RES_UNK,WT_unmapped,RES_unmapped" > benchmarkdna/final_scores2.csv
find benchmarkdna |grep dna.updated_table_with_scores_and_mutations.tsv$ |while read output ; do
    echo $output |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Wildtype |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Wildtype |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Unknown |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Unknown |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    sample=`echo $output|cut -f 2 -d "/" | sed 's/mpf.//'` ;  cat $output |cut -f 3 |grep WT |sort|uniq >$output.wtmappedreadlist ; cat benchmark/$sample.fastq |grep ^"@" |sed 's/^@//' |grep WT >$output.wtallreadlist ; cat $output.wtmappedreadlist $output.wtallreadlist |sort|uniq -c |grep " 1 " |wc -l |tr "\n" ","
    sample=`echo $output|cut -f 2 -d "/" | sed 's/mpf.//'` ;  cat $output |cut -f 3 |grep RES |sort|uniq >$output.resmappedreadlist ; cat benchmark/$sample.fastq |grep ^"@" |sed 's/^@//' |grep RES >$output.resallreadlist ; cat $output.resmappedreadlist $output.resallreadlist |sort|uniq -c |grep " 1 " |wc -l
done >>benchmarkdna/final_scores2.csv

