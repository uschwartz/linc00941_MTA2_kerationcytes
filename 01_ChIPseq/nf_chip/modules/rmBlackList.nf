process rmBL{
  container 'uschwartz/core_docker:v1.0'
  publishDir "${params.outDir}/MacsInFiles/IP", mode: 'copy', pattern: '*.bam'

  input:
  tuple val(sampleID), file(bam), file(bam_idx), file(blacklist)

  output:
  tuple val(sampleID), file("*.bam")

  script:
  """
  bedtools intersect -v -abam $bam -b $blacklist > ${sampleID}_rmBL.bam
  """

}
