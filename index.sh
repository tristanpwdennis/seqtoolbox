for f in `cat list`
do
    samtools index bam/$f".srt.dd.bam"
done
