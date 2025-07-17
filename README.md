# metapointfinder
finds and scores resistance causing mutations in long reads using the information in the AMRFinder database

install:

``git clone https://github.com/aldertzomer/metapointfinder.git``

use: 

``./metapointfinder file.fastq databasefolder outputfolder or ./metapointfinder file.fastq.gz databasefolder outputfolder``

dependencies:

R (msa and Biostrings libs)
wget (to get the AMRFinder databases)
diamond (to align the reads against the reference protein sequences)
kma (to align the reads against the reference DNA sequences
