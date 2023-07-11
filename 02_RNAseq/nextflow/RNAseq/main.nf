#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based RNA-seq pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

 //                           show settings

if (!params.help) {
println ''
 log.info """\

          R N A - s e q   P I P E L I N E
          =============================
          Path variables used in analysis
          reads : ${params.fastqPath}
          output : ${params.outPath}
          STARidx : ${params.STARidxPath}
          GTF.files : ${params.gtfPath}
          Name.Expression: ${params.exprName}


          General options
          pairedEnd : ${params.pairedEnd}
          strandness : ${params.strandness}
          GTFforCounting: ${params.gtfFile}
          Trimming: ${params.trim}

          """.stripIndent()

println ''
}

//                       help message

/*
* Help message
*/
def helpMessage() {
    println ''
    log.info """
    R N A - s e q   P I P E L I N E
    =============================
    Usage:

    nextflow run RNAseq --fastqPath 'path2dirFASTQ' --outPath 'path2output'

    Mandatory arguments:
      --fastqPath           Path to fastq data (must end with .fastq.gz)
      --outPath             Path to write output
      --STARidxPath         Path to STAR index
      --gtfPath             Path to GTF files


    Generic:
      --pairedEnd           Specifiy if input is paired-end sequencing data (allowed: true/false; default: false)
      --strandness          Specifiy directionality of sequencing data (allowed: forward/reverse/unstranded; default: reverse)
      --testRUN             performs a test run of the first 100,000 reads of one sample (default false)
      --gtfFile             Name of GTF file for counting which is located in gtfPath (standard: protein_coding.gtf/all_genes.gtf/protein_coding_and_lincRNA.gtf ;default: protein_coding.gtf)
      --exprName            Expression used to look for fastq files of single end reads (default: '*.fastq.gz')
      --exprNamePE          Expression used to look for fastq files of paired end reads (default: '*{1,2}.fastq.gz')
      --highMemory          Specify if memory overflow in STAR (default: false)


    STAR:
      --zipSTAR             unzipping option used in STAR [default: gunzip]

    Trimming:
      --trim             enable adapter trimming (allowed: true/false; default: false)
      --adapters         Specify the path to the adapter sequence file in fasta format

     """.stripIndent()
     println ''
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// to set strandness for featureCounts and dupRadar
strandness_Map = ["unstranded": 0, "forward" : 1 , "reverse" : 2]

 //                           FASTQC of raw Data

 /*
 * Create a channel for input read files
 */

if(params.testRUN){
        fastq_files=Channel.fromPath("${params.fastqPath}/${params.exprName}")
                .splitFastq(by: 10000, limit: 10000, file: true).first()
} else {
        fastq_files=Channel.fromPath("${params.fastqPath}/${params.exprName}")
}

outPath_defined = (params.testRUN ? "${params.outPath}/testRUN":"${params.outPath}" )

/*
* FASTQC process
*/
process runFastQC{

	publishDir "${outPath_defined}/FastQC_rawData", mode: 'copy'

	input:
	path fastq_input from fastq_files

	output:
	file "*_fastqc.{zip,html}" into FastQC_out

	"""
	fastqc $fastq_input
	"""

}

/*
* multiQC summary
*/
process multiQC{
	publishDir "${outPath_defined}/QC_summaries/00_FastQC_rawData", mode: 'copy'

	input:
	file ('fastqc/*') from FastQC_out.collect().ifEmpty([])

	output:
	file "*multiqc_report.html" into multiqc_report
  file "*_data"

	"""
	multiqc .
	"""
}


 //                           Alignment


 /*
 * get input reads
 */

// optional trimming
if(params.trim){
        if(params.testRUN){
                fastq2trim = ( params.pairedEnd
                                      ? Channel.fromFilePairs("${params.fastqPath}/${params.exprNamePE}", flat:true)
                                      .first().splitFastq(by: 10000_10000, limit: 10000 ,pe:true, file:true)
                                      .map { filename,fwd,rev-> tuple(filename, [fwd,rev]) }
                                      : Channel
                                        .fromPath("${params.fastqPath}/${params.exprName}").first()
                                        .map { file -> tuple(file.simpleName, file.splitFastq(by: 10000, limit: 10000, file: true))}
                                        )

        } else {
                fastq2trim = ( params.pairedEnd
                                      ? Channel.fromFilePairs("${params.fastqPath}/${params.exprNamePE}")
                                      : Channel
                                        .fromPath("${params.fastqPath}/${params.exprName}")
                                        .map { file -> tuple(file.simpleName, file)}
                                        )
        }

        process trim{
                label 'big_mem'
                publishDir "${outPath_defined}/trimming/logs", mode: 'copy', pattern: "*.log"
                input:
                set nameID, file(reads) from fastq2trim
                output:
                set nameID, "*_trimmed.fastq.gz" into fastq2align
                file "*_trimmed.fastq.gz" into fastqcTrim
                file "*.log" into trimLog


                script:
                if (params.pairedEnd) {
                        """
                        java -jar $trimmomatic PE -threads 7 -phred33 $reads ${nameID}_1_trimmed.fastq.gz ${nameID}_1_unpaired.fastq.gz ${nameID}_2_trimmed.fastq.gz ${nameID}_2_unpaired.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 2>${nameID}_trimming.log
                        """
                } else {
                        """
                        java -jar $trimmomatic SE -threads 7 -phred33  $reads ${nameID}_trimmed.fastq.gz  ILLUMINACLIP:${params.adapters}:2:30:10:6 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 2>${nameID}_trimming.log
                        """
                }
        }
        process FastQC_trim{
                publishDir "${outPath_defined}/trimming/FastQC", mode: 'copy'

        	input:
        	path fastq_input from fastqcTrim

        	output:
        	file "*_fastqc.{zip,html}" into FastQC_out_trim

        	"""
        	fastqc $fastq_input
        	"""
        }
        process multiQC_Trim{
        	publishDir "${outPath_defined}/QC_summaries/01_trimming/", mode: 'copy'

        	input:
        	file ('trimmomatic/*') from trimLog.collect().ifEmpty([])
                file ('fastqc/*') from FastQC_out_trim.collect().ifEmpty([])

        	output:
        	file "*multiqc_report.html" into multiqc_report_Trim
                file "*_data"

        	"""
        	multiqc .
        	"""
        }

} else {
        if(params.testRUN){
                fastq2align = ( params.pairedEnd
                                ? Channel.fromFilePairs("${params.fastqPath}/${params.exprNamePE}", flat:true).first()
                                .splitFastq(by: 10000_10000, limit: 10000 ,pe:true, file:true)
                                .map { filename,fwd,rev-> tuple(filename, [fwd,rev]) }
                              : Channel
                                .fromPath("${params.fastqPath}/${params.exprName}").first()
                                .map { file -> tuple(file.simpleName, file.splitFastq(by: 10000, limit: 10000, file: true))}
                                )
        } else {
                fastq2align = ( params.pairedEnd
                              ? Channel.fromFilePairs("${params.fastqPath}/${params.exprNamePE}")
                              : Channel
                                .fromPath("${params.fastqPath}/${params.exprName}")
                                .map { file -> tuple(file.simpleName, file)}
                                )
        }
}





/*
* alignment using STAR
*/
// if error appears in STAR "ERROR: not enough memory for BAM sorting:"

if(params.highMemory){
        process align_unsorted_bam{
          label 'big_mem'

          publishDir "${outPath_defined}/alignment/mapping_reports", mode: 'copy', pattern: "*.final.out"
          input:
          set nameID, file(reads) from fastq2align
          output:
          set nameID, "*.bam" into bam2sort
          file "*.final.out" into STARlog

          script:
          readFilesCmd=(params.testRUN & !params.trim ? '' :"--readFilesCommand $params.zipSTAR -c")
          """
          ### set limit of threads
          ulimit -n 10000

          STAR --outFilterType BySJout \\
                --outFilterMultimapNmax 20 \\
                --alignSJoverhangMin 8 \\
                --alignSJDBoverhangMin 1 \\
                --outFilterMismatchNmax 999 \\
                --alignIntronMin 20 \\
                --alignIntronMax 1000000 \\
                --outFilterMismatchNoverReadLmax 0.04 \\
                --runThreadN 8 \\
                --outSAMtype BAM Unsorted \\
                --outSAMmultNmax 1 \\
                --outMultimapperOrder Random  $readFilesCmd \\
                --genomeDir $params.STARidxPath \\
                --readFilesIn  $reads \\
                --outFileNamePrefix $nameID


          """

        }
        process sortBam{
                label 'big_mem'
                input:
                set nameID, file(reads) from bam2sort
                output:
                set nameID, "*sorted.bam" into bam, bam2bw

                script:
                """
                samtools sort -O BAM -@ 7 -o ${nameID}_sorted.bam $reads
                """
        }

        process bam2bw{
                memory = 70.GB
                publishDir "${outPath_defined}/alignment/profiles/multiple/", mode: 'copy'
                input:
                set nameID, file(reads) from bam2bw
                output:
                file "*.bw"

                script:
                """
                samtools index ${reads}
                bamCoverage --bam ${reads} -o ${nameID}.bw \
                --normalizeUsing CPM
                """
        }

} else {
        process align{
          label 'big_mem'

          publishDir "${outPath_defined}/alignment/mapping_reports", mode: 'copy', pattern: "*.final.out"
          input:
          set nameID, file(reads) from fastq2align
          output:
          set nameID, "*.bam" into bam
          file "*.wig" into wigs
          file "*.final.out" into STARlog

          script:
          readFilesCmd=(params.testRUN & !params.trim ? '' :"--readFilesCommand $params.zipSTAR -c")
          """
          ### set limit of threads
          ulimit -n 10000

          STAR --outFilterType BySJout \\
                --outFilterMultimapNmax 20 \\
                --alignSJoverhangMin 8 \\
                --alignSJDBoverhangMin 1 \\
                --outFilterMismatchNmax 999 \\
                --alignIntronMin 20 \\
                --alignIntronMax 1000000 \\
                --outFilterMismatchNoverReadLmax 0.04 \\
                --runThreadN 8 \\
                --outSAMtype BAM SortedByCoordinate \\
                --outWigType wiggle \\
                --outSAMmultNmax 1 \\
                --outMultimapperOrder Random  $readFilesCmd \\
                --genomeDir $params.STARidxPath \\
                --readFilesIn  $reads \\
                --outFileNamePrefix $nameID


          """

        }
        /*
        * wig to bigwig
        */
        if(!params.testRUN){
                process makeBigWigs{
                        publishDir "${outPath_defined}/alignment/profiles/unique/", mode: 'copy', pattern:"*.Unique.*"
                        publishDir "${outPath_defined}/alignment/profiles/multiple/", mode: 'copy', pattern:"*.UniqueMultiple*"

                        input:
                        file wig from wigs.flatten()

                        output:
                        file "*.bw"

                        script:
                        """
                        wigToBigWig $wig ${params.STARidxPath}/chrNameLength.txt ${wig.baseName}.bw
                        """
                }
        }
}



/*
* multiqc summary STAR alignment
*/
process multiQC_STAR{
	publishDir "${outPath_defined}/QC_summaries/02_alignment", mode: 'copy'

	input:
	file ('STAR/*') from STARlog.collect().ifEmpty([])

	output:
	file "*multiqc_report.html" into multiqc_report_STAR
        file "*_data"

	"""
	multiqc .
	"""
}




//                           Quality control

/*
* Mark read duplicates
*/

process markDuplicates{
        label 'mid_mem'
        errorStrategy 'retry'
        maxRetries 5
        publishDir "${outPath_defined}/alignment/", mode: 'copy', pattern:"*_dupmark.bam*"

        input:
        set nameID, file(bam_file) from bam

        output:
        set nameID, "*_dupmark.bam" into  bamDupMark, bamDupMark_dupRadar, bamCount
        file "*_dupstats.txt" into statsDupMark
        file "*_dupmark.bam.bai"


        script:
        """
        ulimit -c unlimited
        java -Xmx30g -jar /Users/admin/Tools/picard.jar MarkDuplicates \\
         M=${nameID}_dupstats.txt \\
         I=$bam_file \\
         O=${nameID}_dupmark.bam \\
         REMOVE_DUPLICATES=FALSE \\
         OPTICAL_DUPLICATE_PIXEL_DISTANCE=0
         samtools index ${nameID}_dupmark.bam
        """
}

/*
process markDuplicates{
        publishDir "${outPath_defined}/alignment/", mode: 'copy', pattern:"*_dupmark.bam*"

        input:
        set nameID, file(bam_file) from bam

        output:
        set nameID, "*_dupmark.bam" into  bamDupMark, bamDupMark_dupRadar, bamCount
        file "*_dupmark.log" into statsDupMark
        file "*_dupmark.bam.bai"


        script:
        """
        bam dedup --in $bam_file --log ${nameID}_dupmark.log --out ${nameID}_dupmark.bam --noPhoneHome
        samtools index ${nameID}_dupmark.bam
        """
}
*/

/*
* multiqc summary of MarkDuplicates
*/
process multiQC_markDup{
	publishDir "${outPath_defined}/QC_summaries/03_duplicates", mode: 'copy'

	input:
	file ('dupmark/*') from statsDupMark.collect().ifEmpty([])

	output:
	file "*multiqc_report.html" into multiqc_report_dupMark
        file "*_data"

	"""
	multiqc .
	"""
}

/*
* qualimap
*/
process qualimap{
  label 'mid_mem'
        publishDir "${outPath_defined}/alignment/qualimap/", mode: 'copy'

	input:
	set nameID, file(bam_file) from bamDupMark

	output:
        file "${nameID}/*" into qualimap

        script:
        pairedOpt = ( params.pairedEnd ? '-pe':'')
	"""
        export JAVA_HOME=`/usr/libexec/java_home -v 1.8`

	qualimap rnaseq --java-mem-size=16G -a uniquely-mapped-reads -bam $bam_file -gtf ${params.gtfPath}/protein_coding.gtf -outdir ${nameID} $pairedOpt
	"""
}

/*
* multiqc summary of qualimap
*/
process multiQC_qualimap{
	publishDir "${outPath_defined}/QC_summaries/04_qualimap", mode: 'copy'

	input:
	file ('qualimap/*/*') from qualimap.collect().ifEmpty([])

	output:
	file "*multiqc_report.html" into multiqc_report_qualimap
        file "*_data"

	"""
	multiqc .
	"""
}

if(!params.testRUN){

/*
* dupRadar
*/

process dupRadar{
        publishDir "${outPath_defined}/alignment/dupRadar/${nameID}", mode: 'copy'

	input:
	set nameID, file(bam_file) from bamDupMark_dupRadar


	output:
        file "*{pdf,txt}" into dupRadar


        script:
        stranded=strandness_Map[params.strandness]
        pairDup=( params.pairedEnd ? 'TRUE':'FALSE')
	"""
        Rscript ${params.dupRadarPath} $bam_file ${params.gtfPath}/protein_coding.gtf $stranded $pairDup
	"""
}


//                          Count reads
/*
* featureCounts
*/



process countReads{
        publishDir "${outPath_defined}/counts/", mode: 'copy'

	input:
        file bams from bamCount.flatMap { logs, bams -> bams }.collect().ifEmpty([])

	output:
        file "*.txt*" into countOut


        script:
        stranded=strandness_Map[params.strandness]
        pairedOpt = ( params.pairedEnd ? '-p':'')
	"""
        featureCounts $pairedOpt -T 10 -s $stranded  -a ${params.gtfPath}/${params.gtfFile} -o count_table.txt $bams &>count_info.txt
	"""
}

/*
* multiqc summary of counting
*/
process multiQC_counting{
	publishDir "${outPath_defined}/QC_summaries/05_counting", mode: 'copy'

	input:
	file ('featureCounts/*') from countOut.collect().ifEmpty([])

	output:
	file "*multiqc_report.html" into multiqc_report_counting
        file "*_data"

	"""
	multiqc .
	"""
}
}
