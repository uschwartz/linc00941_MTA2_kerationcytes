process fingerprint{
  container 'uschwartz/core_docker:v1.0'
  label 'mid'

  publishDir "${params.outDir}/RUN/04_FingerPrint", mode: 'copy'

  input:
  tuple val(names),file(bams), file(bams_idx), val(group), file(blacklist)

  output:
  file("*_fingerprint.pdf")
  file("*_fingerprintMetrics.txt")

  script:
  //blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  plotFingerprint -b $bams  \
  -plot ${group}_fingerprint.pdf \
  -p $task.cpus \
  --extendReads \
  --outQualityMetrics ${group}_fingerprintMetrics.txt \
  --region "17" \
  --binSize 200 \
  --numberOfSamples 400000 \
  --blackListFileName $blacklist
  """

}

process fingerprint_dedup{
  container 'uschwartz/core_docker:v1.0'
  label 'mid'

  publishDir "${params.outDir}/RUN/05_FingerPrint_rmDUPLICATES", mode: 'copy'

  input:
  tuple val(names),file(bams), file(bams_idx), val(group)

  output:
  file("*_fingerprint_dedup.pdf")
  file("*_fingerprintMetrics_dedup.txt")

  script:
  blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  plotFingerprint -b $bams  \
  -plot ${group}_fingerprint_dedup.pdf \
  -p $task.cpus \
  --extendReads \
  --ignoreDuplicates \
  --outQualityMetrics ${group}_fingerprintMetrics_dedup.txt \
  --region "17" \
  --binSize 200 \
  --numberOfSamples 400000 \
  $blacklistOpt
  """

}
