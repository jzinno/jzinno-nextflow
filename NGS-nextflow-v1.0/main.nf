// Declare syntax version
//nextflow.enable.dsl=2

println "And we're off!..."




ref = file(params.ref)




// Importing Illumina reads form tsv, but the table lacks full paths so we stick them in from a data directory config param
Channel
     .fromPath( params.table )
     .splitCsv(header: ['pair_id', 'read1', 'read2'], sep:'\t')
     .map{ row-> tuple(row.pair_id, file("${params.datadir}" + row.read1), file("${params.datadir}" + row.read2)) }
     .take( 3 )
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
    module load fastp/intel/0.20.1 

    fastp \
    -i ${read1} \
    -I ${read2} \
    -o ${pair_id}_1.fastp.fastq.gz \
    -O ${pair_id}_2.fastp.fastq.gz \
    -h ${pair_id}.fastp.html \
    -j ${pair_id}.fastp.json \
    --length_required 76 \
    --n_base_limit 50 \
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
    module load fastqc/0.11.9

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
    module load bwa/intel/0.7.17
    module load samtools/intel/1.14



    bwa mem \
    -M \
    -t 8 \
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
    module load picard/2.23.8
    module load samtools/intel/1.14



    java -Xmx44g -jar /share/apps/picard/2.23.8/picard.jar SortSam \
        INPUT=${sam_bam[1]} \
        OUTPUT=${pair_id}_sorted.bam \
        SORT_ORDER=coordinate \
        MAX_RECORDS_IN_RAM=10000000 \
        VALIDATION_STRINGENCY=LENIENT

    samtools index ${pair_id}_sorted.bam
    """
}
