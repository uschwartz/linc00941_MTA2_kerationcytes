AnalysisDir=/Volumes/PromisePegasus/_Service_/S026_siChIP_E2F6_MTA2

cd $AnalysisDir


nextflow run ./script/nf_chip/ \
--csvInput 'script/addData/fastqIN.csv' \
--outDir 'run_221117' \
--genomeIdx 'bowtie2idx/bt2_' \
--blacklist $AnalysisDir'/script/addData/blacklisted_regions_comprehensive_ENCFF356LFX_MT_rDNA.bed' \
--genomeDir '/Users/admin/Annotation/GRCh38/bowtie2idx/' \
-resume
