GENOME=~/Genomes/Onchocerca_volvulus/ncbi/oncho_with_wolb.fna

for f in `cat list`
do
    ~/software/gatk/gatk HaplotypeCaller -I bam/$f".srt.dd.bam" -O vcf/$f".g.vcf" -R $GENOME -ERC GVCF
done
