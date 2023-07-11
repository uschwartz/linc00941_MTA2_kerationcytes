AnalysisDir=/Volumes/PromisePegasus/_Service_/S036_skin_MTA_figures

cd $AnalysisDir


nextflow run ./script/Figure4/nf_bigwig/ \
--csvInput 'script/Figure4/addData/bamIN.csv' \
--outDir 'analysis/Figure4/bigwigs'  \
--peaksMTA2 'data/peaks/mta2_filtFDR_narrowPeaks.bed' \
--bamInput 'data/bigwigs/MTA2_bams/pooled_input.bam' \
-resume
