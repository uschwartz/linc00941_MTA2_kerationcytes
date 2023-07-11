process multibigwigsummary_rep{
  container 'uschwartz/core_docker:v1.0'
  label 'big'

  publishDir "${params.outDir}/RUN/07_Rep_on_Peaks/${group}/", mode: 'copy', pattern: '*count.txt'

  input:
  tuple val(group),val(names),file(bw),val(names2),file(bed)
  file(input)

  output:
  tuple val(group),file("*_results.npz")
  file("*_raw_count.txt")

  script:
  """
  multiBigwigSummary BED-file -b $bw $input \
  -o replicates_results.npz \
  --BED $bed \
  --smartLabels \
  --outRawCounts ${group}_raw_count.txt \
  -p $task.cpus
  """

}

process plotCorrelation_rep{
  container 'uschwartz/core_docker:v1.0'

  publishDir "${params.outDir}/RUN/07_Rep_on_Peaks/${group}/", mode: 'copy'

  input:
  tuple val(group),file(corData)

  output:
  file("*_SampleCorrelation.pdf")
  file("*_SampleCorrelation_heat.pdf")
  file("*_SampleCorr_Matrix.tab")

  script:
  """
  plotCorrelation --corData $corData \
  --corMethod pearson \
  --whatToPlot scatterplot \
  -o ${group}_SampleCorrelation.pdf \
  --colorMap RdYlBu \
  --skipZeros \
  --outFileCorMatrix ${group}_SampleCorr_Matrix.tab

  plotCorrelation --corData $corData \
  --corMethod pearson \
  --whatToPlot heatmap \
  -o ${group}_SampleCorrelation_heat.pdf \
  --colorMap RdYlBu \
  --skipZeros \
  """

}

process plotPCA_rep{
  container 'uschwartz/core_docker:v1.0'

  publishDir "${params.outDir}/RUN/07_Rep_on_Peaks/${group}/", mode: 'copy'

  input:
  tuple val(group),file(corData)

  output:
  file("*_SamplePCA.pdf")

  script:
  """
  plotPCA --corData $corData \
  --transpose \
  -o ${group}_SamplePCA.pdf
  """

}


process multibigwigsummary_rep_norm{
  container 'uschwartz/core_docker:v1.0'
  label 'big'

  publishDir "${params.outDir}/RUN/07_Rep_on_Peaks/${group}/", mode: 'copy', pattern: '*count.txt'

  input:
  tuple val(group),val(names), file(bed), val(names2),  file(bw)

  output:
  tuple val(group),file("*_results.npz")
  file("*_norm_count.txt")

  script:
  """
  multiBigwigSummary BED-file -b $bw \
  -o replicates_results.npz \
  --BED $bed \
  --smartLabels \
  --outRawCounts ${group}_norm_count.txt \
  -p $task.cpus
  """

}
