GTF=~/Annotation/GRCh38/Homo_sapiens.GRCh38.106.gtf

OutDir=~/Analysis/S015_skin_RNAseq/NFrun/counts
mkdir $OutDir"/newerAnno"

bamDir=~/Analysis/S015_skin_RNAseq/NFrun/alignment

cd $bamDir

featureCounts -T 10 -s 2  -a $GTF \
 -o $OutDir"/newerAnno"/count_table.txt *.bam &>$OutDir"/newerAnno"/count_info.txt
