process fastqc{
  container "quay.io/biocontainers/fastqc:0.11.9--0"
  publishDir "${params.outDir}/QC/01_FASTQC/${sampleID}", mode: 'copy', pattern: "*.html"

  input:
  tuple val(sampleID), file(read)

  output:
  file "*_fastqc.zip"
  file "*_fastqc.html"

  script:
  """
  fastqc $read
  """
}
