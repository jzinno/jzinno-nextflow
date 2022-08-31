//Begin Workflow
println "And we're off!..."



//importing refrence sequence from params, this is not really needed beacuse you can just reference params.ref but I kept it
ref = file(params.ref)




// Importing Illumina reads form tsv, but the table lacks full paths so we stick them in from a data directory config param
Channel
     .fromPath( params.table )
     .splitCsv(header: ['pair_id', 'read1', 'read2'], sep:'\t')
     .map{ row-> tuple(row.pair_id, file("${params.datadir}" + row.read1), file("${params.datadir}" + row.read2)) }
     .set { read_pairs_ch }




//run fastp 
process fastp {
    publishDir "${params.out}/trimmed_reads", mode:'copy'

    input:
    set pair_id, file(read1), file(read2) from read_pairs_ch

    output:
    set val(pair_id), file("${pair_id}_{1,2}.fastp.fastq.gz") into trimmed_reads_ch 
    set val(pair_id), file("${pair_id}.fastp.html"), file("${pair_id}.fastp.json") into qc_ch
    
    script:
    """
    #!/bin/bash 
    module purge
    

    /gpfs/home/jpz239/fastp \
    -i ${read1} \
    -I ${read2} \
    -o ${pair_id}_1.fastp.fastq.gz \
    -O ${pair_id}_2.fastp.fastq.gz \
    -h ${pair_id}.fastp.html \
    -j ${pair_id}.fastp.json \
    --length_required 76 \
    --n_base_limit 50 \
    --dedup \
    --detect_adapter_for_pe
    """
}

//Splitting up this channel for alignment and fastqc
trimmed_reads_ch.into { fastqc_in; alignment_in }

//run fastqc
process fastqc {
    publishDir "${params.out}/QC", mode:'copy'
    
    input:
    set pair_id, file(reads) from fastqc_in

    output:
    set val(pair_id), file("${pair_id}_{1,2}.fastp_fastqc.html") into fastqc_reports_ch 


    script:
    """
    #!/bin/bash
    module purge
    module load fastqc/0.11.7

    fastqc ${reads[0]} ${reads[1]}
    """
}

//run bwa mem
process bwamem {
    publishDir "${params.out}/aligned", mode:'copy'

    input:
    set pair_id, file(reads) from alignment_in

    output:
    set val(pair_id), file("${pair_id}.{sam,bam}") into mappings_ch
    script:
    """
    #!/bin/bash
    module purge
    module load bwa/0.7.17
    module load samtools/1.14



    bwa mem \
    -M \
    -t 8 \
    -R "@RG\\tID:${pair_id}\\tPL:ILLUMINA\\tSM:${pair_id}" \
    ${ref} \
    ${reads[0]} \
    ${reads[1]} > ${pair_id}.sam

    samtools view -b -h ${pair_id}.sam > ${pair_id}.bam
    """
}


//run Picard coordiante sorting
process picard {
    publishDir "${params.out}/sorted", mode:'copy'

    input:
    set pair_id, file(sam_bam) from mappings_ch

    output:
    set val(pair_id), file("${pair_id}_sorted.bam"), file("${pair_id}_sorted.bam.bai") into sorted_ch

    script:
    """
    #!/bin/bash
    module purge
    module load picard/2.18.11
    module load samtools/1.14



    java -Xmx44g -jar /gpfs/share/apps/picard/2.18.11/libs/picard.jar SortSam \
        INPUT=${sam_bam[1]} \
        OUTPUT=${pair_id}_sorted.bam \
        SORT_ORDER=coordinate \
        MAX_RECORDS_IN_RAM=10000000 \
        VALIDATION_STRINGENCY=LENIENT

    samtools index ${pair_id}_sorted.bam
    """
}


//run bqsr report
process bqsr_report {
    publishDir "${params.out}/bqsr_reports", mode:'copy'

    input:
    set pair_id, file(bam) from sorted_ch

    output:
    set val(pair_id), file("${bam}"), file("${pair_id}_bqsr.report") into bqsr_reports_ch

    script:
    """
    #!/bin/bash
    module purge

    module load gatk/4.2.1.0


    gatk --java-options "-Xmx4g" BaseRecalibrator  \
        -I ${bam} \
        -R ${ref} \
        --known-sites /gpfs/data/kirchhofflab/Resources/hg38/gatk_bundle/dbsnp_146.hg38.vcf.gz \
        --known-sites /gpfs/data/kirchhofflab/Resources/hg38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -O ${pair_id}_bqsr.report

    """
}

//run bqsr apply
process bqsr_apply {
    publishDir "${params.out}/bqsr_bams", mode:'copy'

    input:
    set pair_id, file(bam), file(report) from bqsr_reports_ch

    output:
    set val(pair_id), file("${pair_id}_bqsr.bam"), file("${pair_id}_bqsr.bam.bai") into prepped_bams_ch

    script:
    """
    #!/bin/bash
    module purge

    module load gatk/4.2.1.0
    module load samtools/1.14

    gatk --java-options "-Xmx4g" ApplyBQSR \
        -I ${bam} \
        -R ${ref} \
        --bqsr-recal-file ${report} \
        -O ${pair_id}_bqsr.bam


    samtools index ${pair_id}_bqsr.bam
    
    """
}

//run haplotype_caller
process  haplotype_caller {
    publishDir "${params.out}/gVCFs", mode:'copy'

    input:
    set pair_id, file(bam) from prepped_bams_ch

    output:
    set val(pair_id), file("${pair_id}.g.vcf.gz") into gvcf_ch

    script:
    """
    #!/bin/bash
    module purge

    module load gatk/4.2.1.0
    
    gatk --java-options "-Xmx10g" HaplotypeCaller  \
      -R ${ref} \
      -I ${bam} \
      -O ${pair_id}.g.vcf.gz \
      -ERC GVCF

    """
}