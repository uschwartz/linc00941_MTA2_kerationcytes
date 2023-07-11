AnalysisDir=~/Analysis/S015_skin_RNAseq/

cd $AnalysisDir

nextflow run ~/00_scripts/nextflow/RNAseq/main.nf  \
	--fastqPath $AnalysisDir/data/pooled \
	--outPath $AnalysisDir/NFrun \
	--STARidxPath ~/Annotation/GRCh38/STARidx \
	--gtfPath ~/Annotation/GRCh38/nextflow \
	--gtfFile all_genes.gtf
