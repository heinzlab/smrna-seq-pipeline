#!/usr/bin/env nextflow

/*
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quality control
 - 2:   BBDuk for adapter trimming
 - 3.1: Bowtie 1 alignment against miRBase mature miRNA
 - 3.2: Post-alignment processing of miRBase mature miRNA counts
 - 3.3: edgeR analysis on miRBase mature miRNA counts
        - TMM normalization and a table of top expression mature miRNA
        - MDS plot clustering samples
        - Heatmap of sample similarities
 - 4.1: Bowtie 1 alignment against miRBase hairpin for the unaligned reads in step 3
 - 4.2: Post-alignment processing of miRBase hairpin counts
 - 4.3: edgeR analysis on miRBase hairpin counts
        - TMM normalization and a table of top expression hairpin
        - MDS plot clustering samples
        - Heatmap of sample similarities
 - 5.1: Bowtie 2 alignment against host reference genome
 - 5.2: Post-alignment processing of Bowtie 2
 - 6:   NGI-Visualization of Bowtie 2 alignment statistics
 - 7:   MultiQC
----------------------------------------------------------------------------------------
*/

def helpMessage() {
	log.info"""
	==============================================
	smRNA-seq Pipeline
	==============================================
	Usage:

	nextflow run smrna-seq-pipeline/main.nf --reads '*.fastq.gz' --genome hg38

	Mandatory arguments:
	 --reads                       Path to input data (must be surrounded with quotes).
	 --genome                      Name of iGenomes reference.

	Trimming options:
	--length [int]                Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18

	Other options:
	--outdir                      The output directory where the results will be saved. Default: results
	--skip_qc                     Skip all QC steps aside from MultiQC.
	--skip_fastqc                 Skip FastQC.

	References:
	--saveReference               Save the generated reference files the the Results directory.
	""".stripIndent()
}

// Show help message
params.help = false
if (params.help){
	helpMessage()
	exit 0
}

// Reference path configuration
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bt2index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.bt2indices = null
params.mature = params.genome ? params.genomes[ params.genome ].mature ?: false : false
params.hairpin = params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
params.adapters = params.genome ? params.genomes[ params.genome ].adapters ?: false: false

// Validate inputs
if( !params.mature || !params.hairpin ){
	exit 1, "Missing mature / hairpin reference indexes! Is --genome specified?"
}
if( params.adapters ){
	adapters = file(params.adapters)
	if ( !adapters.exists() ) exit 1, "adapters file not found: ${params.adapters}"
}
if( params.mature ){
	mature = file(params.mature)
	if( !mature.exists() ) exit 1, "Mature file not found: ${params.mature}"
}
if( params.hairpin ){
	hairpin = file(params.hairpin)
	if( !hairpin.exists() ) exit 1, "Hairpin file not found: ${params.hairpin}"
}
if( params.gtf ){
	gtf = file(params.gtf)
	if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}
if( params.bt2index ){
	bt2_index = file("${params.bt2index}.fa")
	bt2_indices = Channel.fromPath( "${params.bt2index}*.bt2" ).toList()
	if( !bt2_index.exists() ) exit 1, "Reference genome Bowtie 2 not found: ${params.bt2index}"
} else if( params.bt2indices ){
	bt2_indices = Channel.from(params.readPaths).map{ file(it) }.toList()
}

// Create a channel for input reads
if(params.readPaths){
    Channel
        .from(params.readPaths)
        .map { file(it) }
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .into { raw_reads_fastqc; raw_reads_bbduk }
} else {
    Channel
        .fromPath( params.reads )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
        .into { raw_reads_fastqc; raw_reads_bbduk }
}

// Header log info
log.info """=======================================================
smRNA-seq pipeline - svenner lab
======================================================="""
def summary = [:]
summary['Reads']               = params.reads
summary['Genome']              = params.genome
summary['Trim Min Length']     = params.length
summary['miRBase mature']      = params.mature
summary['miRBase hairpin']     = params.hairpin
if(params.bt2index)            summary['Bowtie2 Index'] = params.bt2index
if(params.gtf)                 summary['GTF Annotation'] = params.gtf
summary['adapters']            = params.adapters
summary['Output dir']          = params.outdir
summary['Working dir']         = workflow.workDir
summary['Container']           = params.container
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Script dir']          = workflow.projectDir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================="

/*
 * PREPROCESSING - Build Bowtie index for mature and hairpin
 */
process makeBowtieIndex {

    publishDir path: { params.saveReference ? "${params.outdir}/bowtie/reference" : params.outdir },
               saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file mature from mature
    file hairpin from hairpin

    output:
    file 'mature_idx.*' into mature_index
    file 'hairpin_idx.*' into hairpin_index

    script:
    """
    fasta_formatter -w 0 -i $mature -o mature_igenome.fa
    fasta_nucleotide_changer -d -i mature_igenome.fa -o mature_idx.fa
    bowtie-build mature_idx.fa mature_idx
    fasta_formatter -w 0 -i $hairpin -o hairpin_igenome.fa
    fasta_nucleotide_changer -d -i hairpin_igenome.fa -o hairpin_idx.fa
    bowtie-build hairpin_idx.fa hairpin_idx
    """
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$reads"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    file reads from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - BBDuk
 */
process bbduk {
    tag "$reads"
    publishDir "${params.outdir}/bbduk", mode: 'copy'

    input:
    file reads from raw_reads_bbduk
    file adapters from adapters

    output:
    file '*.gz' into trimmed_reads_bowtie, trimmed_reads_bowtie2, trimmed_reads_insertsize

    script:
    tg_length = "--length ${params.length}"
    c_r1 = params.clip_R1 > 0 ? "--clip_R1 ${params.clip_R1}" : ''
    tpc_r1 = params.three_prime_clip_R1 > 0 ? "--three_prime_clip_R1 ${params.three_prime_clip_R1}" : ''
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bbduk.sh in=$reads out=${prefix}_trimmed.fastq.gz ref=$adapters -ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20
    """
}

/*
 * STEP 2.1 - Insertsize
 */

process insertsize {
    tag "$reads"
    publishDir "${params.outdir}/bbduk/insertsize", mode: 'copy'

    input:
    file reads from trimmed_reads_insertsize

    output:
    file '*.insertsize' into insertsize_results

    script:
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(zcat $reads) >${prefix}.insertsize
    """
}


/*
 * STEP 3 - Bowtie against miRBase mature miRNA
 */
process bowtie_miRBase_mature {
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_mature", mode: 'copy', pattern: '*.mature_unmapped.fq.gz'

    input:
    file reads from trimmed_reads_bowtie
    file index from mature_index

    output:
    file '*.mature.bam' into miRBase_mature_bam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bowtie \\
        $index_base \\
        -q <(zcat $reads) \\
        -p 2 \\
        -t \\
        -k 1 \\
        -m 1 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un ${prefix}.mature_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.mature.bam
    gzip ${prefix}.mature_unmapped.fq
    """
}

/*
 * STEP 4 - Bowtie against miRBase hairpin
 */
process bowtie_miRBase_hairpin {
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_hairpin", mode: 'copy', pattern: '*.hairpin_unmapped.fq.gz'

    input:
    file reads from mature_unmapped_reads
    file index from hairpin_index

    output:
    file '*.hairpin.bam' into miRBase_hairpin_bam
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - '.mature_unmapped.fq.gz'
    """
    bowtie \\
        $index_base \\
        -p 2 \\
        -t \\
        -k 1 \\
        -m 1 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        -q <(zcat $reads) \\
        --un ${prefix}.hairpin_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.hairpin.bam
    gzip ${prefix}.hairpin_unmapped.fq
    """
}


/*
 * STEP 5 - Post-alignment processing for miRBase mature and hairpin
 */
def wrap_mature_and_hairpin = { file ->
    if ( file.contains("mature") ) return "miRBase_mature/$file"
    if ( file.contains("hairpin") ) return "miRBase_hairpin/$file"
}

process miRBasePostAlignment {
    tag "$input"
    publishDir "${params.outdir}/bowtie", mode: 'copy', saveAs: wrap_mature_and_hairpin

    input:
    file input from miRBase_mature_bam.mix(miRBase_hairpin_bam)

    output:
    file "${input.baseName}.count" into miRBase_counts
    file "${input.baseName}.sorted.bam" into miRBase_bam
    file "${input.baseName}.sorted.bam.bai" into miRBase_bai

    script:
    """
    samtools sort ${input.baseName}.bam -o ${input.baseName}.sorted.bam
    samtools index ${input.baseName}.sorted.bam
    samtools idxstats ${input.baseName}.sorted.bam > ${input.baseName}.count
    """
}

/*
 * STEP 6 - edgeR miRBase feature counts processing
 */
process edgeR_miRBase {
    publishDir "${params.outdir}/edgeR", mode: 'copy', saveAs: wrap_mature_and_hairpin

    input:
    file input_files from miRBase_counts.toSortedList()

    output:
    file '*.{txt,pdf}' into edgeR_miRBase_results

    script:
    """
    edgeR_miRBase.r $input_files
    """
}

/*
 * STEP 7.1 and 7.2 IF A GENOME SPECIFIED ONLY!
 */
if( params.gtf && params.bt2index) {

    /*
     * STEP 7.1 - Bowtie 2 against reference genome
     */
    process bowtie2 {
        tag "$reads"
        publishDir "${params.outdir}/bowtie2", mode: 'copy'

        input:
        file reads from trimmed_reads_bowtie2
        file bt2_indices

        output:
        file '*.bowtie2.bam' into bowtie2_bam, bowtie2_bam_for_unmapped

        script:
        index_base = bt2_indices[0].toString()  - ~/\.\d+\.bt2/
        prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        bowtie2 \\
            -x $index_base \\
            -U $reads \\
            -k 10 \\
            --very-sensitive \\
            -p 8 \\
            -t \\
            | samtools view -bT $index_base - > ${prefix}.bowtie2.bam
        """
    }

    /*
     * STEP 7.2 - Bowtie 2 Statistics about unmapped reads against ref genome
     */

    process bowtie2_unmapped {
        tag "${input_files[0].baseName}"
        publishDir "${params.outdir}/bowtie2/unmapped", mode: 'copy'

        input:
        file input_files from bowtie2_bam_for_unmapped.toSortedList()

        output:
        file 'unmapped_refgenome.txt' into bowtie2_unmapped

        script:
        """
        for i in $input_files
        do
          printf "\${i}\t"
          samtools view -c -f0x4 \${i}
        done > unmapped_refgenome.txt
        """
    }


    /*
     * STEP 7.3 - NGI-Visualizations of Bowtie 2 alignment statistics
     */
    process ngi_visualizations {
        tag "$bowtie2_bam"
        publishDir "${params.outdir}/bowtie2/ngi_visualizations", mode: 'copy'

        input:
        file gtf from gtf
        file bowtie2_bam

        output:
        file '*.{png,pdf}' into bowtie2_ngi_visualizations

        script:
        // Note! ngi_visualizations needs to be installed!
        // See https://github.com/NationalGenomicsInfrastructure/ngi_visualizations
        """
        #!/usr/bin/env python
        from ngi_visualizations.biotypes import count_biotypes
        count_biotypes.main('$gtf','$bowtie2_bam')
        """
    }

}
