process fdr_peak_filter{
  container 'uschwartz/core_docker:v1.0'
  publishDir "${params.outDir}/RUN/06_CallPeaks/${group}/", mode: 'copy'

  input:
  tuple  val(group), file(bed), val(names)

  output:
  tuple  val(names),file("*_narrowPeaks.bed"), val(group)

  script:
  """
  awk -v OFS='\t' '\$9 >= 20 {print}' $bed >>${group}_filtFDR_narrowPeaks.bed
  """

}
