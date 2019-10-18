mkdir -p bam/bamsummarystats

for f in `cat list`

do
    samtools flagstat bam/$f".srt.dd.bam" > bam/bamsummarystats/$f".flagfile"

    genomeCoverageBed -ibam bam/$f".srt.dd.bam" > bam/bamsummarystats/$f".flagfile"

done
