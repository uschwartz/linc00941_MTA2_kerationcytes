process computeMatrix{
  container 'uschwartz/core_docker:v1.0'
  label 'mid'

  input:
  tuple val(group),val(names), file(peak), val(names2),  file(bw)

  output:
  tuple val(group), file("*_compMx.gz")

  script:
  """
  computeMatrix reference-point \
  -S $bw \
  -R $peak \
  -o  ${group}_compMx.gz\
  -p $task.cpus \
  --referencePoint center \
  -b 1000 -a 1000 \
  --smartLabels
  """
}

process plotHeatmap{
        container 'uschwartz/core_docker:v1.0'
        publishDir "${params.outDir}/RUN/07_Rep_on_Peaks/${group}/", mode: 'copy'

        input:
        tuple val(group), file(compmx)

        output:
        file("*.pdf")

        script:
        """
        plotHeatmap  -m $compmx \
        -o replicates_heatmap.pdf \
        -y "logFC_over_input"
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
