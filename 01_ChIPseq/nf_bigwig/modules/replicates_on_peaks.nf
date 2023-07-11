process computeMatrix{
  container 'uschwartz/core_docker:v2.0'
  label 'mid'

  input:
  file(bw)
  file(peak)

  output:
  file("*_compMx.gz")

  script:
  """
  computeMatrix reference-point \
  -S $bw \
  -R $peak \
  -o  rep_compMx.gz\
  -p $task.cpus \
  --referencePoint center \
  -b 1500 -a 1500 \
  --smartLabels
  """
}

process plotHeatmap{
        container 'uschwartz/core_docker:v2.0'
        publishDir "${params.outDir}/Rep_on_MTA2_Peaks/", mode: 'copy'

        input:
        file(compmx)

        output:
        tuple file("*.pdf"), file("rawMatrixHeat.gz")

        script:
        """
        plotHeatmap  -m $compmx \
        -o replicates_heatmap.pdf \
        --colorMap Blues \
        -y "norm. MTA2 ChIP coverage" \
        -x "distance to MTA2 peak center" \
        --refPointLabel "MTA2 peak center" \
        --outFileNameMatrix rawMatrixHeat.gz \
        --heatmapWidth 6 --heatmapHeight 10 \
        --legendLocation "none"
        """
}


process plotHeatmap_cluster{
        container 'uschwartz/core_docker:v1.0'
        publishDir "${params.outDir}/RUN/07_Rep_on_Peaks/${group}/", mode: 'copy'

        input:
        tuple val(group), file(compmx)

        output:
        file("*.pdf")

        script:
        """
        plotHeatmap  -m $compmx \
        --kmeans 3 \
        -o replicates_heatmap_cluster.pdf \
        -y "logFC_over_input"
        """
}
