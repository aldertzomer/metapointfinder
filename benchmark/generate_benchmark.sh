# this generates a set of benchmark fastq files based on the amrprot database of amrfinder as processed by metapointfinder and then processes these and scores them for TP, TN, FP
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
	echo "wgsim -N 1000000 -1 100 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.100.$r.fastq /dev/null"
	echo "wgsim -N 500000 -1 200 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.200.$r.fastq /dev/null"
	echo "wgsim -N 100000 -1 1000 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.1000.$r.fastq /dev/null"
	echo "wgsim -N 20000 -1 5000 -d 0 -S 42 -e 0 -r $r benchmark/$cut.6000.fasta benchmark/$cut.5000.$r.fastq /dev/null"
    done
done |parallel -j 64 >/dev/null

# run metapointfinder on all outputs. change the locations of the tool and the db to suit your needs
for cut in Eecoli Estaau Emyctu ; do
    for r in 0.00 0.01 0.02 0.03 0.05 0.1 0.2 ; do
	echo "/mnt/data/tools/metapointfinder/metapointfinder/metapointfinder benchmark/$cut.100.$r.fastq /mnt/data/db/metapointfinderdb benchmark/mpf.$cut.100.$r"
	echo "/mnt/data/tools/metapointfinder/metapointfinder/metapointfinder benchmark/$cut.200.$r.fastq /mnt/data/db/metapointfinderdb benchmark/mpf.$cut.200.$r"
	echo "/mnt/data/tools/metapointfinder/metapointfinder/metapointfinder benchmark/$cut.1000.$r.fastq /mnt/data/db/metapointfinderdb benchmark/mpf.$cut.1000.$r"
	echo "/mnt/data/tools/metapointfinder/metapointfinder/metapointfinder benchmark/$cut.5000.$r.fastq /mnt/data/db/metapointfinderdb benchmark/mpf.$cut.5000.$r"
    done
done |parallel -j 64

# this thing checks if the information in the header is the same as the prediction. this counts unique hits. We only check WT reads and Wildtype prediction (TN) and RES reads and Resistance prediction (TP). There will be many false negatives and unknowns as many of the shorter reads do not cover the resistance determining region, therefore these are not scored. False Positives can occur by chance in the reads with mutation.
echo "file,TP,TN,FP" > benchmark/final_scores.csv
find benchmark |grep prot.updated_table_with_scores_and_mutations.tsv |while read output ; do
    echo $output |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Wildtype |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}'|sort|uniq |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l
done >>benchmark/final_scores.csv

# this thing checks if the information in the header is the same as the prediction. This counts all hits. We only check WT reads and Wildtype prediction (TN) and RES reads and Resistance prediction (TP). There will be many false negatives and unknowns as many of the shorter reads do not cover the resistance determining region, therefore these are not scored. False Positives can occur by chance in the reads with mutation.
echo "file,TP,TN,FP" > benchmark/final_scores2.csv
find benchmark |grep prot.updated_table_with_scores_and_mutations.tsv |while read output ; do
    echo $output |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Wildtype |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_RES ; done |wc -l |tr "\n" ","
    cat $output  |sed 's/_0:/\t/' |cut -f 1,3,11 |awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,"_");if(length(a)>2){$2=a[1];for(i=2;i<=length(a)-2;i++)$2=$2"_"a[i]}print}' |grep Resistant |while read drug name result ; do  echo $name |grep "$drug"_WT ; done |wc -l
done >>benchmark/final_scores2.csv