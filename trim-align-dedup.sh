#####
#Specify paths
#####

#To BAM fi
BAMDIR=~/onchogenome/faust_schisto/sequencing-data-pipeline/bam/
FASTQDIR=~/onchogenome/faust_schisto/sequencing-data-pipeline/fastq/
VCFDIR=~/onchogenome/faust_schisto/sequencing-data-pipeline/bam/vcf/
GENOME=~/Genomes/schisto/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa
TRIMMODIR=~/software/Trimmomatic-0.39/
THREADS=6

#####
#Prepare and index genome
#####

echo "bwa indexing genome...."

#bwa index $GENOME

echo "samtools indexing genome..."

#samtools faidx $GENOME

#start loop over sample list

echo "starting loop...."

for f in `cat list`

do

######
#trim reads
######

#java -jar $TRIMMODIR/trimmomatic-0.39.jar PE $FASTQDIR/untrim/$f"_R1_001.fastq.gz" $FASTQDIR/untrim/$f"_R2_001.fastq.gz" $FASTQDIR/trim/$f"_R1.trim.fq.gz" $FASTQDIR/trim/$f"_R1.trim.unpaired.fq.gz" $FASTQDIR/trim/$f"_R2.trim.fq.gz" $FASTQDIR/trim/$f"_R2.trim.unpaired.fq.gz" ILLUMINACLIP:$TRIMMODIR/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

######
#align reads
######

bwa mem -t $THREADS $GENOME $FASTQDIR/trim/$f"_R1.trim.fq.gz" $FASTQDIR/trim/$f"_R2.trim.fq.gz" | samtools sort > $BAMDIR/$f".srt.bam"

######
#remove duplicates
######

java -jar ~/software/picard/build/libs/picard.jar MarkDuplicates I=$BAMDIR/$f".srt.bam" O=$BAMDIR/$f".srt.dd.bam" M=$BAMDIR/$f".M"

done
