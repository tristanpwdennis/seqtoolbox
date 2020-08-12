# seqtoolbox
sequencing utilities for DNDI-BMGF funded project on Onchocerca volvulus

## Base

1. run init.sh
2. download fastq files to fastq/untrim
3. fiddle with the file extensions and paths in main script (will need to update paths based on where all the software is)
4. make sure bwa mem and samtools are in your PATH and gatk/picard .jars are compiled
5. run

## Utility scripts
1. vcftofasta.sh run as `bash vcftofasta.sh <vcf.vcf.gz>`
2. will take specified vdf, tabix index and convert to a multifasta file

## Analysis scripts
R scripts for:
-Visualising uniformity and depth of sequencing coverage
-Modelling the effect of coverage on kinship assignment, and plotting kinship network

1. Open in RStudio
2. Tinker with paths
3. Run

Python environment and notebook for faffing about with sequencing data
1. Install and activate environment with python >= 3.5
2. Clone repo
3. `cd path/to/repo`
4. `pip3 install -r requirements.txt`
5. `pip3 install -e .`
6. `jupyter notebook`
7. Import vcf using first bit and then alter according to what kind of parsing you want to do
