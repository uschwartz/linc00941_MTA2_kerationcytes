process alignment{
  container 'uschwartz/core_docker:v1.0'
  label 'big'
  publishDir "${params.outDir}/RUN/01_ALIGNMENT", mode: 'copy', pattern: '*_aligned.{bam,bam.bai}'
  publishDir "${params.outDir}/QC/02_ALIGNMENT/", mode: 'copy', pattern: '*_alignment_stats.txt'


  input:
  tuple val(sampleID), file(read1), file(read2), path(bw2idx)

  output:
  file "*_alignment_stats.txt"
  tuple val(sampleID), file("*_aligned.bam"), file("*_aligned.bam.bai")

  script:
  """
  bowtie2 -t \
  --threads $task.cpus \
  --very-sensitive-local \
  --no-discordant \
  --no-mix \
  --dovetail \
  -x $params.genomeIdx \
  -1 $read1 \
  -2 $read2 \
  2> ${sampleID}_alignment_stats.txt \
  | samtools view -bS -q 30 -f 2 -@ $task.cpus - | samtools sort -@ $task.cpus - > ${sampleID}"_aligned.bam"

  samtools index -b ${sampleID}"_aligned.bam"
  """

}
