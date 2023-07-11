process bam2bw{
  container 'uschwartz/core_docker:v2.0'
  label 'mid'
  publishDir "${params.outDir}/bw_scaleFactor", mode: 'copy'


  input:
  tuple val(sampleID), file(bam), val(sizeFac)


  output:
  tuple val(sampleID), file("*.bw")

  script:
  """
  samtools index $bam

  bamCoverage -b $bam \
   -o ${sampleID}".bw" \
   -p $task.cpus \
   --scaleFactor $sizeFac \
   --extendReads \
  """

}

process bam2bw_pooled{
  container 'uschwartz/core_docker:v2.0'
  label 'mid'
  publishDir "${params.outDir}/bw_pooled", mode: 'copy'


  input:
  file(bam)
  val(name)


  output:
  file("*.bw")

  script:
  """
  samtools index $bam

  bamCoverage -b $bam \
   -o ${name}"_pooled.bw" \
   -p $task.cpus \
   --extendReads \
   --normalizeUsing 'CPM'
  """

}
