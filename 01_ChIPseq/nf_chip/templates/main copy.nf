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

 //                           show settings
 if (!params.help) {
         include{settings} from './modules/setting'
         settings()
 }

 /*/                       help message
 // Show help message
 if (params.help) {
     include{helpMessage} from './modules/help'
     helpMessage()
     exit 0
 }
*/


 //                      workflow
// read csv file
// forward reads
Channel.fromPath(params.csvInput)
        .splitCsv(header:true)
        .map{ row -> tuple(row.name,file(row.path_fwdReads))}
        .set{samples_fwd_ch}

// reverse reads
Channel.fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.name,file(row.path_revReads))}
      .set{samples_rev_ch}


//Channel for fastqc
samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}
//Channel for alignment
samples_fwd_ch.join(samples_rev_ch).set{samplePair_ch}

//////////// grouping on antibody ///////////////
Channel.fromPath(params.csvInput)
       .splitCsv(header:true)
       .map{ row -> tuple(row.name,row.antibody, row.input)}
       .set{groups}
// get IP and corresponding Input
//IPs
groups.map{
 name,antibody,input -> tuple(antibody.toLowerCase(),name)
}.filter { it[0]!='input'}.groupTuple(by:0).set{ips}
// Input names
groups.map{
 name,antibody,input -> tuple(antibody.toLowerCase(),input)
}.filter { it[0]!='input'}.groupTuple(by:0).set{inputs}
//sample and group ID
ips.join(inputs).map{
 name,ip,input -> tuple(ip.plus(input).unique(), name)
 }.transpose().set{ip_groups}


 groups.map{
   name,antibody,input -> tuple(name,antibody.toLowerCase())
 }.filter { it[1]=='input'}
 .map{ name,antibody -> name }.set{input_ids}

 groups.map{
    name,antibody,input -> tuple(name,antibody.toLowerCase())
  }.filter { it[1]!='input'}.set{ips_ids}


////////////////////////////////////////////////////////////


// load modules
//fastqc
include{fastqc} from './modules/raw_qc'
//multiqc
include{multiqc} from './modules/multiqc'
//alignment to ref genome
include{alignment} from './modules/align'
//visualize bam
include{bam2bw; bam2bw_dedup} from './modules/bam2bw'
//correlate samples
include{multibigwigsummary; plotCorrelation; plotPCA} from './modules/correlate_samples'
include{multibigwigsummary_dedup; plotCorrelation_dedup; plotPCA_dedup} from './modules/correlate_samples'
//fingerprint
include{fingerprint_dedup; fingerprint} from './modules/fingerprint'

//remove blacklisted regions
include{rmBL} from './modules/rmBlackList'
//pool inputs
include{pool} from './modules/pool_input'
// macs callpeak
include{macs} from './modules/peakcalling'
// input normalized
include{norm_input} from './modules/norm_input_bw'
// filter peaks
include{fdr_peak_filter} from './modules/filterPeaks'
// replicates
include{computeMatrix; plotHeatmap} from './modules/replicates_on_peaks'
include{multibigwigsummary_rep; plotCorrelation_rep; plotPCA_rep; multibigwigsummary_rep_norm} from './modules/correlate_replicates'

//refine peak peakcalling
include{refinePeaks;refinePeaks_merge} from './modules/workflow_refinePeakCalling'

workflow{
  fastqc(sampleSingle_ch)
  alignment(samplePair_ch)
  bam2bw(alignment.out[1])

  //correlate
  multibigwigsummary(bam2bw.out.map{name,bw -> file(bw)}.collect())
  plotCorrelation(multibigwigsummary.out)
  plotPCA(multibigwigsummary.out)
  //fingerprint
  fingerprint(alignment.out[1].join(ip_groups).groupTuple(by:3))

  // remove rmDuplicates
  bam2bw_dedup(alignment.out[1])
  multibigwigsummary_dedup(bam2bw_dedup.out.map{name,bw -> file(bw)}.collect())
  plotCorrelation_dedup(multibigwigsummary_dedup.out)
  plotPCA_dedup(multibigwigsummary_dedup.out)
  //fingerprint
  fingerprint_dedup(alignment.out[1].join(ip_groups).groupTuple(by:3))

  //rmBlackList
  rmBL(alignment.out[1])
  //pool Input
  pool(rmBL.out.join(input_ids).map{names,bam -> bam}.collect())

  // Visualize input normed
  norm_input{rmBL.out.join(ips_ids).combine(pool.out)}

  //call peaks
  macs(rmBL.out.join(ips_ids).groupTuple(by:2).combine(pool.out))

  //filter peaks
  fdr_peak_filter(macs.out[1])

  //check reproducibility
  norm_input.out.groupTuple(by:2).set{group_normbw}
  computeMatrix(fdr_peak_filter.out.join(group_normbw,by:2))
  plotHeatmap(computeMatrix.out[0])

  bam2bw.out.join(input_ids).map{ name,bam-> bam}.toList().set{bw_inputs}
  multibigwigsummary_rep(bam2bw.out.join(ips_ids).groupTuple(by:2).join(fdr_peak_filter.out,by:2),bw_inputs)
  plotCorrelation_rep(multibigwigsummary_rep.out[0])
  plotPCA_rep(multibigwigsummary_rep.out[0])

  multibigwigsummary_rep_norm(fdr_peak_filter.out.join(group_normbw,by:2))

  multiqc(fastqc.out[0].mix(alignment.out[0]).collect())

  //----------------- refinePeaks----------------
  // select samples
  Channel.of(['E2F4_1', 'e2f4_v2'], ['E2F4_2', 'e2f4_v2'],
        ['E2F6_2', 'e2f6_v2'], ['E2F6_3', 'e2f6_v2']).set{ips_select}

  refinePeaks(rmBL.out.join(ips_select).groupTuple(by:2).combine(pool.out),
  group_normbw, bam2bw.out, bw_inputs)

  // select samples E2F6 individual peak calling
  Channel.of(['E2F6_2', 'e2f6_v3_rep2'], ['E2F6_3', 'e2f6_v3_rep3']).set{ips_select2}
  refinePeaks_merge(rmBL.out.join(ips_select2).groupTuple(by:2).combine(pool.out),
  group_normbw, bam2bw.out, bw_inputs)
  
}
