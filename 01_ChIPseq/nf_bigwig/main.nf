#!/usr/bin/env nextflow

 /*
 ===============================================================================
             nextflow based ChIP-seq explore analysis pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2



 //                      workflow


////////////// read csv file  ////////////
// forward reads
Channel.fromPath(params.csvInput)
        .splitCsv(header:true)
        .map{ row -> tuple(row.name,file(row.path_bam),row.sizeFactor)}
        .set{samples_bam}

////////////// read peak file  ////////////
Channel.fromPath(params.peaksMTA2).set{mta2_peaks}
////////////// read peak file  ////////////
Channel.fromPath(params.bamInput).set{input_bam}

 ////////////////////////////////////////////////////////////

include {bam2bw;bam2bw_pooled} from './modules/bam2bw'
// replicates at MTA2 peaks
include{computeMatrix; plotHeatmap} from './modules/replicates_on_peaks'
//correlate samples
include{multibigwigsummary; plotCorrelation; plotPCA} from './modules/correlate_samples'
//pool
include{pool} from './modules/pool'
//
include{bw_log2FC} from './modules/bw_log2FC'


workflow{
        bam2bw(samples_bam)

        //replicates at MTA2 peaks

        bam2bw.out.map{names,bw -> names}.toSortedList().flatten()
        .join(bam2bw.out).map{names,bw -> file(bw)}.collect().set{sort_bw}

        computeMatrix(sort_bw,mta2_peaks)
        plotHeatmap(computeMatrix.out)


        //check reproducibility
        //correlate
        multibigwigsummary(bam2bw.out.map{name,bw -> file(bw)}.toSortedList())
        plotCorrelation(multibigwigsummary.out)
        plotPCA(multibigwigsummary.out)

        // pool all replicates to MTA2
        pool(samples_bam.map{
                name,bam,sizeFac -> file(bam)
                }.collect(),'MTA2')
        bw_log2FC(pool.out,input_bam,'MTA2')
        bam2bw_pooled(pool.out,'MTA2')


}
