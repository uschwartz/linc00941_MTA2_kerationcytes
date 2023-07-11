// macs callpeak
include{macs} from './peakcalling'
// filter peaks
include{fdr_peak_filter} from './filterPeaks'
// replicates
include{computeMatrix; plotHeatmap; plotHeatmap_cluster} from './replicates_on_peaks'
include{multibigwigsummary_rep; plotCorrelation_rep} from './correlate_replicates'


workflow refinePeaks{
        take:
        input_macs
        group_normbw
        bam2bw
        bw_inputs

        main:
        //call peaks
        macs(input_macs)
        //filter peaks
        fdr_peak_filter(macs.out[1])

        // get sample selection
        group_normbw.transpose().join(fdr_peak_filter.out.transpose())
        .groupTuple(by:3)
        .map{
          names,bw,group_v1,peak,group_v2 -> tuple(group_v2[0],names,peak,group_v1,bw)
          }.set{compMX}

        computeMatrix(compMX)
        plotHeatmap(computeMatrix.out[0])
        plotHeatmap_cluster(computeMatrix.out[0])

        bam2bw.join(fdr_peak_filter.out.transpose()).groupTuple(by:3)
        .map{
          names,bw,bed,group -> tuple(group,names,bw,names,bed[0])
          }.set{corr}

        multibigwigsummary_rep(corr,bw_inputs)
        plotCorrelation_rep(multibigwigsummary_rep.out[0])
}

process mergePeaks{
  container 'uschwartz/core_docker:v1.0'
  publishDir "${params.outDir}/RUN/06_CallPeaks/${group}/", mode: 'copy'

  input:
  tuple val(names), file(bed1), file(bed2), val(group)

  output:
  tuple  val(names),file("*_peaks.bed"), val(group)

  script:
  """
  bedtools intersect -a $bed1 -b $bed2 >intersect_peaks.bed
  """
}


workflow refinePeaks_merge{
        take:
        input_macs
        group_normbw
        bam2bw
        bw_inputs

        main:
        //call peaks
        macs(input_macs)
        //filter peaks
        fdr_peak_filter(macs.out[1])

        macs.out[1].toList().map{
          it1,it2 -> tuple(tuple(it1[2],it2[2]).flatten(),it1[1],it2[1], "e2f6_v3_intersect")
          }.set{intersect_peaks}

        // intersect
        mergePeaks(intersect_peaks)

        // get sample selection


        group_normbw.transpose().join(mergePeaks.out.transpose())
        .groupTuple(by:3)
        .map{
          names,bw,group_v1,peak,group_v2 -> tuple(group_v2[0],names,peak,group_v1,bw)
          }.set{compMX}

        computeMatrix(compMX)
        plotHeatmap(computeMatrix.out[0])
        plotHeatmap_cluster(computeMatrix.out[0])

}
