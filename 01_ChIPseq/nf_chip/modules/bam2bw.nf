process bam2bw{
  container 'uschwartz/core_docker:v1.0'
  label 'mid'
  publishDir "${params.outDir}/RUN/02_VisualizeBAM/all", mode: 'copy'


  input:
  tuple val(sampleID), file(bam), file(bam_idx), file(blacklist)
 

  output:
  tuple val(sampleID), file("*.bw")

  script:
  //blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  bamCoverage -b $bam \
   -o ${sampleID}".bw" \
   -p $task.cpus \
   --normalizeUsing 'CPM' \
   --extendReads \
   --blackListFileName $blacklist
  """

}


process bam2bw_dedup{
  container 'uschwartz/core_docker:v1.0'
  label 'mid'
  publishDir "${params.outDir}/RUN/02_VisualizeBAM/rmDuplicates/", mode: 'copy'


  input:
  tuple val(sampleID), file(bam), file(bam_idx), file(blacklist)


  output:
  tuple val(sampleID), file("*.bw")

  script:
  //blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  bamCoverage -b $bam \
   -o ${sampleID}"_dedup.bw" \
   -p $task.cpus \
   --normalizeUsing 'CPM' \
   --extendReads \
   --ignoreDuplicates \
   --blackListFileName $blacklist
  """

}
