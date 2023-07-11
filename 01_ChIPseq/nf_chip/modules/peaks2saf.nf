process peaks2saf{
  container 'uschwartz/core_docker:v1.0'
  publishDir "${params.outDir}/RUN/06_CallPeaks/${group}/", mode: 'copy'


  input:
  tuple  val(names),file(peak), val(group)


  output:
  tuple  val(group),file("*_filt_peaks.saf")

  script:
  """
  echo -e "GeneID\tChr\tStart\tEnd\tStrand" >${group}"_filt_peaks.saf"
  awk -v OFS='\t' '{print \$4,\$1,\$2,\$3,"."}' \
   $peak >> ${group}"_filt_peaks.saf"
  """

}
