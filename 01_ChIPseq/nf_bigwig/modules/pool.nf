process pool{
  label 'big'
  container 'uschwartz/core_docker:v2.0'


  input:
  file(bams)
  val(name)
  
  output:
  file("*.bam")


  script:
  """
  samtools merge -@ $task.cpus pooled_${name}.bam $bams
  """
}
