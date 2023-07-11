process bw_log2FC{
  container 'uschwartz/core_docker:v2.0'
  label 'mid'
  publishDir "${params.outDir}//bw_log2FC/", mode: 'copy'


  input:
  file(bam)
  file(input)
  val(name)

  output:
  tuple val(name), file("*.bw")

  script:
  """
  samtools index -@ $task.cpus -b $bam
  samtools index -@ $task.cpus -b $input

  bamCompare -b1 $bam -b2 $input \
  -o ${name}_log2FC.bw \
  --extendReads \
  --pseudocount 10 \
  -p $task.cpus
  """

}
