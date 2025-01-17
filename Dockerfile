FROM buchfink/diamond
# check diamond version used for metapointfinder

RUN  apt-get install -y wget curl

RUN wget -c https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa
RUN wget -c https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt-mutation.tsv
#replace slashes by underscores
RUN cat AMRProt-mutation.tsv | sed 's/\//_/g' > AMRProt-mutation_underscore.tsv
#remove this line and replace class with a file with only the AB class you want to check
RUN cat AMRProt-mutation_underscore.tsv | cut -f 6 | sort | uniq > class

RUN cat class |while read class ; do \
    cat AMRProt-mutation_underscore.tsv | cut -f 2,4,6 | grep -w $class | sort -nk 1 | while read accession mutation class ; do \
	if [ "$accession" != "$previousaccession" ] ; then \
	    echo $class $previousaccession $changes_str | sed 's/,$//' \
	    changes_str=`` \
	fi \
	mutation_cleaned=`echo $mutation | cut -f 2 -d _` \
	changes_str=`echo "$mutation_cleaned"",""$changes_str"` \
	previousaccession=`echo $accession` \
    done \
done |tr " " "\t" > AMRProt-mutation_underscore_combined.tsv

RUN diamond makedb --in AMRProt.fa --db AMRProt 1>diamondmakedb.log 2>diamondmakedb.error

FROM rocker/r-ver

RUN R -e "install.packages('msa')"
RUN R -e "install.packages('Biostrings')"

