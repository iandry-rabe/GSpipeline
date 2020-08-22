#!/usr/bin/env nextflow

// v.1.0.0

// General architecture of a nextflow run
// /data/analysis_name
// ├── Input
// |   └── Raw FastQ folder
// ├── References
// │   ├── GATK References
// │   ├── bwa References
// │   └── Panels
// │       ├── Panels csv data
// │       ├── Panels bed (pharmaco, covid)
// │       └── Reports templates
// ├── Logs
// ├── Output
// │   ├── FastQ
// │   ├── FastQC
// │   ├── BAM
// │   ├── VCF
// │   └── Panels
// │       ├── csv
// │       ├── ng ??
// │       └── pdf
// ├── Scripts
// ├── Temporary
// └── Tools

// analysis name and output prefixes
params.run_name = null

// input folder
params.input_folder = null

// results folder
params.results_folder = null

// parameters
params.cpus = null
params.genome_bwa = null
params.genome_gatk = null
params.genome_basename = null
params.gatk_dbsnp = null
params.gatk_indels = null


// Nextflow script path
nfPipelinePath = "${workflow.scriptFile}"

// input
analysis_name = params.run_name

// outputs
//analysis_folder = file(params.results_folder + "/analysis")
analysis_folder = file("analysis")
analysis_folder.mkdir()

//analysis_folder = file(params.results_folder + "/fastq")
fastq_folder = file("fastq")
fastq_folder.mkdir()

// check fastq files and store into channels (R1 and R2)
listOfFiles = file(params.input_folder)
R1Fq = Channel.from(listOfFiles).filter(~/.*(_1.fq|R1.fastq|_R1_(\d+).fastq).gz$/)
R2Fq = Channel.from(listOfFiles).filter(~/.*(_2.fq|R2.fastq|_R2_(\d+).fastq).gz$/)

genome_bwa = Channel.fromPath(params.genome_bwa + '*')

// cat fastq
process catFq1 {
	publishDir "${fastq_folder}/", pattern: "*.gz"
	publishDir "${fastq_folder}/", pattern: "*.ok"

	input: 
	file R1Fq

	output:
	file("${analysis_name}.R1.fastq.gz") into R1FqSample
	file("catR1.ok") into catR1_ok

	script:
	"""
	echo $R1Fq
	cat $R1Fq > ${analysis_name}.R1.fastq.gz

	echo "cat R1 ok" > catR1.ok
	"""
}

// cat fastq
process catFq2 {
	publishDir "${fastq_folder}/", pattern: "*.gz"
	publishDir "${fastq_folder}/", pattern: "*.ok"

	input: 
	file R2Fq
	file catR1_ok

	output:
	file("${analysis_name}.R2.fastq.gz") into R2FqSample
	file("catR2.ok") into catR2_ok

	script:
	"""
	echo $R2Fq
	cat $R2Fq > ${analysis_name}.R2.fastq.gz

	echo "cat R2 ok" > catR2.ok
	"""
}

// launch bwa
process bwaMapping {
	publishDir "${analysis_folder}", pattern: "*.bam"

	input:
	file R1FqSample
	file R2FqSample
	file catR2_ok
	file genome_bwa

	output:
	file("*.bam") into bam_files

	script:
	"""
	if [ -s ${R2FqSample} ]; then
		/pipeline/tools/bwa/bwa mem -t ${params.cpus} ${params.genome_basename} ${R1FqSample} ${R2FqSample} | \\
		/pipeline/tools/samtools/samtools view -bS -q1 -@ ${params.cpus} > ${analysis_name}.bam
	else
		/pipeline/tools/bwa/bwa mem -t ${params.cpus} ${params.genome_basename} ${R1FqSample} | \\
		/pipeline/tools/samtools/samtools view -bS -q1 -@ ${params.cpus} > ${analysis_name}.bam
	fi
	"""
}

// launch gatk

process addReadGroups {
	publishDir "${analysis_folder}", pattern: "*.bam"
	publishDir "${analysis_folder}", pattern: "*.bai"

	input:
	file bam_files

	output:
	file("*_RG.bam") into bam_RG
	file("*_RG.bam.bai") into bam_RG_index

	script:
	"""
	java -jar build/lib/picardcloud.jar AddOrReplaceReadGroups \\
	I=${bam_files}} \\
	O=${analysis_name}_RG.bam \\
	RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
	"""
}

process markDups {
	publishDir "${analysis_folder}", pattern: "*.bam"
	publishDir "${analysis_folder}", pattern: "*.bai"
	publishDir "${analysis_folder}", pattern: "*.dups"

	input:
	file bam_RG

	output:
	file("*_RG.dedup.bam") into (bam_dedup_1, bam_dedup_2)
	file("*_RG.dedup.bam.bai") into bam_dedup_index
	file("*_markDuplicate.dups") into bam_duplicates_dups

	script:
	"""
	gatk MarkDuplicatesSpark \\
	-I ${bam_RG} \\
	-O ${analysis_name}_RG.dedup.bam \\
	-M ${analysis_name}_markDuplicate.dups \\
	--remove-all-duplicates \\
	--create-output-bam-index
	"""
}

process baseRecal {
	publishDir "${analysis_folder}", pattern: "*.bam"
	publishDir "${analysis_folder}", pattern: "*.bai"
	publishDir "${analysis_folder}", pattern: "*.table"

	input:
	file bam_dedup_1

	output:
	file("recalibration.table") into recal_table

	script:
	"""
	gatk BaseRecalibrator \\
	-I ${bam_dedup} \\
	-R ${params.genome_gatk} \\
	--known-sites ${params.gatk_dbsnp} \\
	--known-sites ${params.gatk_indels} \\
	-O recalibration.table
	"""
}

process applyBQSR {
	publishDir "${analysis_folder}", pattern: "*.bam"
	publishDir "${analysis_folder}", pattern: "*.bai"

	input:
	file bam_dedup_2
	file recal_table

	output:
	file("*_RG.dedup.recal.bam") into bam_recal
	file("*_RG.dedup.recal.bam.bai") into bam_recal_index

	script:
	"""
	gatk ApplyBQSR \\
	-R ${params.genome_gatk} \\
	-I ${bam_dedup_2} \\
	--bqsr-recal-file ${recal_table} \\
	-O ${analysis_name}_RG.dedup.recal.bam \\
	--create-output-bam-index
	"""
}

process haplotypeCaller {
	publishDir "${analysis_folder}", pattern: "*.vcf"
	publishDir "${analysis_folder}", pattern: "*.vcf.idx"

	input:
	file bam_recal

	output:
	file("*.vcf") into vcf_all
	file("*.vcf.idx") into vcf_all_index

	script:
	"""
	gatk --java-options \"-Xmx8g\" HaplotypeCaller \\
	-R ${params.genome_gatk}} \\
	-I ${bam_recal} \\
	-O ${analysis_name}.all.vcf \\
	-D ${params.gatk_dbsnp}
	"""
}



///////// add createdir