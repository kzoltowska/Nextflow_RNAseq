#!/usr/bin/env nextflow
//Pipeline to process raw sequencing data from a public repository (GEO)
//To run the pipeline: get the SRR numbers of the samples you want to process
//Get the .gtf and .fna files of the respective organism
// General parameters
params.datadir = "${launchDir}/data"
params.outdir = "${launchDir}/results"

//SRR IDs to use for download
//Include in the "" the SRR numbers that you want to process
list_SRR = ["", ""]

//Fastq-dump parameters
//either "" for SE or --split3 for PE
params.SE_PE = "" 
//how many reads to download, for testing purposes only it is 10000
//remove when want to download entire file
params.X = "10000"

//Trimmomatic parameters
params.threads=4
params.slidingwindow="SLIDINGWINDOW:4:15"
params.leadtrain="LEADING:3 TRAILING:3"
params.minlen="MINLEN:36"
params.adapters="ILLUMINACLIP:TruSeq3-SE:2:30:1"
params.average_quality="AVGQUAL:30"

//STAR parameters
params.threadsindex = 12
params.genomeSAindexNbases = 10
params.sjdbOverhang = 36
params.outSAMtype = "BAM SortedByCoordinate"

//HTSEQ_count parameters
params.strandness = "reverse"
params.type = "exon"
params.mode = "union"
params.format = "bam"
params.minq = 10
params.idattr = "gene_id"

srr_ch = channel.fromList(list_SRR)

process fastqdump {
publishDir "${params.datadir}", mode: 'copy', overwrite: true
   
   input:
    val SRR

   output:  
    path("SRR*.fastq.gz"), emit: reads
    path("versions.txt")

    script:
    """
    echo `fastq-dump --version` > versions.txt
    fastq-dump --gzip ${SRR} -X ${params.X}
    
    """
}

process fastqc {
      publishDir "${params.outdir}/fastqc_raw", mode: 'copy', overwrite: true
      container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
   
   input:
    path(reads)
    
   output:  
    path("*{html,zip}"), emit: fastqc_raw
    path("versions.txt")

    script:
    """
    echo `fastqc --version` > "versions.txt"
    fastqc ${reads}
    """
}

process trimmomatic {
      publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true
      container 'quay.io/biocontainers/trimmomatic:0.35--6'

      input:
      path(reads)

      output:
      path("${reads}1_trimmed.fastq.gz"), emit: trimmed
      path("versions.txt")

      script:
      """
      echo `trimmomatic --version` > "versions.txt"
      trimmomatic SE -phred33 ${reads} ${reads}1_trimmed.fastq.gz \\
      -threads ${params.threads} ${params.slidingwindow} ${params.leadtrain} ${params.minlen} \\
      ${params.average_quality}
      """
}

process fastqc_trimmed {
      publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy', overwrite: true
      container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
   input:
    path(reads_t)

   output:  
    path("*.{html,zip}"), emit: fastqc_trimmed
    path("versions.txt")

    script:
    """
      echo `fastqc --version` > "versions.txt"
    fastqc ${reads_t}
    """
}

process multiqc_raw {
      publishDir "${params.outdir}/fastqc_raw", mode: 'copy', overwrite: true
      container 'quay.io/biocontainers/multiqc:1.9--py_1'
   
   input:
    path(inputfiles)

   output:  
    path("*.html"), emit: multiqc
    path("versions_m.txt")

    script:
    """
    echo `multiqc --version` > "versions_m.txt"
    multiqc ${inputfiles}
    """
}

process multiqc_trimmed {
      publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy', overwrite: true
      container 'quay.io/biocontainers/multiqc:1.9--py_1'
   
   input:
    path(inputfiles_t)

   output:  
    path("*.html"), emit: multiqc_trimmed
    path("versions_m.txt")

    script:
    """
    echo `multiqc --version` > "versions_m.txt"
    multiqc ${inputfiles_t}
    """
}

//These channels use previously downloaded files placed here in the genomes directory
//.gtf and .fna
//in the example pipeline for the speed of execution yeast files were used
gtf_ch=Channel.value("/home/kasia/nextflow/Nextflow_RNASeq/genomes/genomicSC.gtf")
genome_ch=Channel.value("/home/kasia/nextflow/Nextflow_RNASeq/genomes/genomeSC.fna")

process star_idx {
    publishDir "${params.outdir}/index_dir"
    container "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0"

    input:
    path genome
    path gtf
    
    output:
    path "index_dir/", emit: index
    path("versions.txt")

    script:
    """
    echo `STAR --version` > "versions.txt"
    mkdir index_dir
    
    STAR --runThreadN ${params.threadsindex} \\
      --runMode genomeGenerate \\
      --genomeDir index_dir/ \\
      --genomeFastaFiles ${genome} \\
      --genomeSAindexNbases ${params.genomeSAindexNbases} \\
      --sjdbGTFfile ${gtf} \\
      --sjdbOverhang ${params.sjdbOverhang} 
    """
}

process star_alignment {
     publishDir "${params.outdir}/mapped-reads/", mode: 'copy', overwrite: true   
     label 'high'
     container "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0"

    input:
     path(reads_to_map) 
     path indexDir
     path gtf

     output:
      path("*.bam"), emit: align_bam
      path("*ReadsPerGene.out.tab"), emit: star_count
      path("versions.txt")

       script:

    """
    echo `STAR --version` > "versions.txt"
    STAR  \\
        --readFilesIn <(zcat ${reads_to_map})\\
        --runThreadN ${params.threadsindex} \\
        --outSAMtype ${params.outSAMtype} \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix ${reads_to_map} \\
        --genomeDir ${indexDir} \\
        --quantMode GeneCounts
    """
}

process htseq {
      publishDir "${params.outdir}/counts/", mode: 'copy', overwrite: true   
     label 'high'
     container "quay.io/biocontainers/htseq:2.0.5--py38h8c35140_1"

      input:
      path(mapped_reads)
      path gtf

      output:
      path("*.txt"), emit: counts
      path("versions.txt")

      script:
      """
      echo `htseq-count --version` > "versions.txt"
      htseq-count --stranded ${params.strandness} --type ${params.type} \\
      --mode ${params.mode} --format ${params.format} -a ${params.minq} --idattr ${params.idattr} \\
      ${mapped_reads} ${gtf} > "${mapped_reads}.txt"
      """
      }

workflow {
      fastqdump(srr_ch)
      fastqc(fastqdump.out.reads)
      trimmomatic(fastqdump.out.reads)
      fastqc_trimmed(trimmomatic.out.trimmed)
      multiqc_raw(fastqc.out.fastqc_raw.collect())
      multiqc_trimmed(fastqc_trimmed.out.fastqc_trimmed.collect())
      star_idx(genome_ch, gtf_ch)
      star_alignment(trimmomatic.out.trimmed, star_idx.out.index.collect(), gtf_ch)
      htseq(star_alignment.out.align_bam, gtf_ch)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Time to complete workflow execution: $workflow.duration"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: $workflow.errorMessage"
}