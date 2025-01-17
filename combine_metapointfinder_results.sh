#!/bin/bash
# usage bash combine_metapointfinder_results.sh regexp , e.g. bash combine_metapointfinder_results.sh GyrA or bash combine_metapointfinder_results.sh QUINOLONE
echo sample class WT R
ls *.updated_table_with_scores_and_mutations.tsv |sed 's/.updated_table_with_scores_and_mutations.tsv//' |while read sample ; do
    echo $1 | while read class ; do
	WT=`cat $sample.updated_table_with_scores_and_mutations.tsv |grep -v class |grep $class|cut -f 7 |grep -w 0 |wc -l`
	R=`cat $sample.updated_table_with_scores_and_mutations.tsv |grep -v class |grep $class|cut -f 7 |grep -v -w 0 |wc -l`
	echo $sample $class $WT $R
    done
done