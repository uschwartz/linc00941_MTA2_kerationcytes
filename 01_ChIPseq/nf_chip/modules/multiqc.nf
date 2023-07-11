process multiqc{
  container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
  publishDir "${params.outDir}/QC/multiqc/", mode: 'copy'

  input:
  file('*')
  output:
  file "*multiqc_report.html"
  file "*_data"

  script:
  """
  multiqc .
  """
}
