process norm_input{
  container 'uschwartz/core_docker:v1.0'
  label 'mid'
  publishDir "${params.outDir}/RUN/02_VisualizeBAM/inputNorm/", mode: 'copy'


  input:
  tuple val(name),file(bam), val(group), file(input)

  output:
  tuple val(name), file("*.bw"), val(group) 

  script:
  """
  samtools index -@ $task.cpus -b $bam
  samtools index -@ $task.cpus -b $input

  bamCompare -b1 $bam -b2 $input \
  -o ${name}_logFC.bw \
  --extendReads \
  -p $task.cpus
  """

}
