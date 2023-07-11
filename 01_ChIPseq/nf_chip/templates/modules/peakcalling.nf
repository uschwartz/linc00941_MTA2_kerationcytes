process macs{
  container 'quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4'
  publishDir "${params.outDir}/RUN/06_CallPeaks/${group}/", mode: 'copy'


  input:
  tuple val(names),file(bams),  val(group), file(input)

  output:
  file("*")
  tuple val(group), file("*.narrowPeak"), val(names)

  script:
  """
  macs2 callpeak \
  -t $bams \
  -c $input \
  -n $group \
  -f BAMPE \
  -g 2.7e9 \
  --keep-dup auto \
  --cutoff-analysis
  """

}
