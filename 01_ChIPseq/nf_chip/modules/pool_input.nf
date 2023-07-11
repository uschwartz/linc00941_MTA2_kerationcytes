process pool{
  label 'big'
  container 'uschwartz/core_docker:v1.0'
  publishDir "${params.outDir}/MacsInFiles/Input", mode: 'copy'

  input:
  file(bams)
  output:
  file("pooled_input.bam")


  script:
  """
  samtools merge -@ $task.cpus pooled_input.bam $bams
  """
}
