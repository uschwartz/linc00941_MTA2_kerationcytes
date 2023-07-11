process subsample_bam{
  container 'uschwartz/core_docker:v1.0'
  label 'big'

  input:
  tuple  val(name),file(bam),file(idx),val(group), val(sub)

  output:
  tuple  val(group),val(name),file("*.bam"),val(sub)

  script:
  """
  samtools view -s $sub -b  -@ $task.cpus $bam >$name"_"$sub".bam"
  """

}


process count_tags_sub{
  container 'uschwartz/core_docker:v1.0'
  label 'big'
  publishDir "${params.outDir}/RUN/09_subsample/sub_${sub}/${group}/", mode: 'copy'


  input:
  tuple  val(group),val(name),file(bams),val(sub),file(saf)

  output:
  tuple val(group),val(sub),file("${group}_${sub}.csv")


  script:
  """
  featureCounts  -F SAF -a $saf \
        -o ${group}_${sub}.csv \
        --fracOverlap 0.8 \
        -T $task.cpus -p --largestOverlap \
        $bams
  """

}

process sigNumb_sub{
  container 'uschwartz/r_core:v2.0'
  publishDir "${params.outDir}/RUN/09_subsample/sub_${sub}/${group}/", mode: 'copy'


  input:
  tuple val(group),val(sub),file(csv)

  output:
  file("*summary.txt")
  file("*.rda")


  script:
  """
  get_signif_genes.R
  """

}
