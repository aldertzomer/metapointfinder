# metapointfinder
finds and scores resistance causing mutations in long reads using the information in the AMRFinder database

install:

``git clone https://github.com/aldertzomer/metapointfinder.git``

use: 

``./metapointfinder file.fastq or ./metapointfinder file.fastq.gz``

dependencies:
R (msa and Biostrings libs)
wget (to get the AMRFinder databases)
diamond (to align the reads against the reference sequences)


