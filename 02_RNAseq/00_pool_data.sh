
AnalysisDir=~/Analysis/S015_skin_RNAseq/
cd $AnalysisDir"/data"

mkdir pooled

for entry in raw/*L005*.fastq.gz; do

	echo $entry
	preID=$(echo $entry | cut -d'/' -f 2-)
	mainID=$(echo $preID | cut -d'_' -f 1-4)

	cat "raw/"$mainID*.fastq.gz > pooled/$mainID".fastq.gz"
done
