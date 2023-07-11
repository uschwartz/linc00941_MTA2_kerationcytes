process multibigwigsummary{
  container 'uschwartz/core_docker:v2.0'
  label 'big'

  input:
  file(bw)

  output:
  file("*_results.npz")

  script:
  """
  multiBigwigSummary bins -b $bw  \
  -o all_results.npz \
  --smartLabels \
  -bs 500 \
  -p $task.cpus
  """

}

process plotCorrelation{
  container 'uschwartz/core_docker:v2.0'

  publishDir "${params.outDir}/SampleCorrelation", mode: 'copy'

  input:
  file(corData)

  output:
  file("*_SampleCorrelation.pdf")
  file("*_SampleCorr_Matrix.tab")

  script:
  """
  plotCorrelation --corData $corData \
  --corMethod spearman \
  --whatToPlot heatmap \
  -o all_SampleCorrelation.pdf \
  --colorMap RdBu \
  --skipZeros \
  --outFileCorMatrix all_SampleCorr_Matrix.tab
  """

}

process plotPCA{
  container 'uschwartz/core_docker:v2.0'

  publishDir "${params.outDir}/SampleCorrelation", mode: 'copy'

  input:
  file(corData)

  output:
  file("*_SamplePCA.pdf")
  file("replicates_PCA.tab")

  script:
  """
  plotPCA --corData $corData \
  --transpose \
  -o all_SamplePCA.pdf \
  --outFileNameData replicates_PCA.tab
  """

}


process multibigwigsummary_dedup{
  container 'uschwartz/core_docker:v1.0'
  label 'big'

  input:
  file(bw)

  output:
  file("*_results.npz")

  script:
  """
  multiBigwigSummary bins -b $bw  \
  -o rmDups_results.npz \
  --smartLabels \
  -p $task.cpus
  """

}

process plotCorrelation_dedup{
  container 'uschwartz/core_docker:v1.0'

  publishDir "${params.outDir}/RUN/03_SampleCorrelation", mode: 'copy'

  input:
  file(corData)

  output:
  file("*_SampleCorrelation.pdf")
  file("*_SampleCorr_Matrix.tab")

  script:
  """
  plotCorrelation --corData $corData \
  --corMethod spearman \
  --whatToPlot heatmap \
  -o rmDups_SampleCorrelation.pdf \
  --colorMap RdYlBu \
  --skipZeros \
  --outFileCorMatrix rmDups_SampleCorr_Matrix.tab
  """

}

process plotPCA_dedup{
  container 'uschwartz/core_docker:v1.0'

  publishDir "${params.outDir}/RUN/03_SampleCorrelation", mode: 'copy'

  input:
  file(corData)

  output:
  file("*_SamplePCA.pdf")

  script:
  """
  plotPCA --corData $corData \
  --transpose \
  -o rmDups_SamplePCA.pdf
  """

}
