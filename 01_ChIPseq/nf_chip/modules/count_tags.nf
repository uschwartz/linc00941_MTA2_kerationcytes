process count_tags{
  container 'uschwartz/core_docker:v1.0'
  label 'big'
  publishDir "${params.outDir}/RUN/08_Count/${group}/", mode: 'copy'


  input:
  tuple  val(group),val(names),file(bams),file(saf)

  output:
  file("*_readCounts.csv")
  file("*readCounts.csv.summary")


  script:
  """
  featureCounts  -F SAF -a $saf \
        -o ${group}_readCounts.csv \
        --fracOverlap 0.8 \
        -T $task.cpus -p --largestOverlap \
        $bams
  """

}
