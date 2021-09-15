#!/usr/bin/env nextflow

/**
===============================
Exome Pipeline
===============================

This Pipeline runs VEP on a set of VCF files

### Homepage / git
git@github.com:ikmb/vep.git
### Implementation
Re-Implemented in Q3 2021

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Pipeline version

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
IKMB VEP pipeline | version ${params.version}
===============================================================================
Usage: nextflow run ikmb/vep --assembly GRCh38 --vcfs /path/to/*.vcf.gz

Required parameters:
--assembly                     Name of the reference assembly to use
--email 		       Email address to send reports to (enclosed in '')
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

// #############
// INPUT OPTIONS
// #############

// Giving this pipeline run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

if (params.run_name == false) {
	log.info "No run name was specified, using ${run_name} instead"
}

// This will eventually enable switching between multiple assembly versions
// Currently, only hg19 has all the required reference files available
params.assembly = "GRCh38"

FASTA = file(params.genomes[ params.assembly ].fasta)
DBSNP = file(params.genomes[ params.assembly ].dbsnp )

Channel.fromPath(params.vcfs)
	.ifEmpty { exit 1; "No VCF files found" }
	.map { b -> [ file(b), file("${b}.tbi") ] }
	.set { vcfs }

// Whether to send a notification upon workflow completion
params.email = false

if (params.assembly != "GRCh38") {
	exit 1, "VEP is not currently set up to work with defunct assembly versions...please use GRCh38"
}

if (!params.vep_cache_dir || !params.vep_plugin_dir) {
	exit 1, "Missing VEP cache and/or plugin directory..."
}

if (params.dbnsfp_db) {
	dbNSFP_DB = file(params.dbnsfp_db)
	if (!dbNSFP_DB.exists()) {
		exit 1, "Could not find the specified dbNSFP database..."
	}
}

if (params.dbscsnv_db) {
	dbscSNV_DB = file(params.dbscsnv_db)
	if ( !dbscSNV_DB.exists() ) {
		exit 1, "Could not find the specified dbscSNV database..."
	}
} else {
	exit 1, "No dbscSNV database defined for this execution profile..."
}

if (params.cadd_snps  && params.cadd_indels ) {
	CADD_SNPS = file(params.cadd_snps)
	CADD_INDELS = file(params.cadd_indels)
	if (!CADD_SNPS.exists() || !CADD_INDELS.exists() ) {
		exit 1, "Missing CADD SNPs and/or Indel references..."
	}
} else {
	exit 1, "CADD SNP and/or Indel reference files not defined for this execution profile..."
}

// Header log info
log.info "========================================="
log.info "VEP pipeline v${params.version}"
log.info "Nextflow Version:             $workflow.nextflow.version"
log.info "Assembly version:             ${params.assembly}"
log.info "-----------------------------------------"
log.info "Command Line:                 $workflow.commandLine"
log.info "Run name:                     ${run_name}"
if (workflow.containerEngine) {
        log.info "Container engine:             ${workflow.containerEngine}"
}
log.info "========================================="


// WORKFLOW starts here

process vep {

        label 'vep'

	publishDir "${params.outdir}/VEP", mode: 'copy'

	input:
        set file(vcf),file(vcf_index) from vcfs

	output:
        file(vcf_vep)

	script:
	vcf_vep = vcf.getBaseName() + ".vep.vcf.gz"

	"""
                vep --offline \
                --cache \
                --dir ${params.vep_cache_dir} \
                --species homo_sapiens \
                --assembly $params.assembly \
                -i $vcf \
                --format vcf \
                -o $vcf_vep --dir_plugins ${params.vep_plugin_dir} \
                --plugin dbNSFP,$dbNSFP_DB,${params.dbnsfp_fields} \
                --plugin dbscSNV,$dbscSNV_DB \
                --plugin CADD,${params.cadd_snps},${params.cadd_indels} \
                --plugin ExACpLI \
		--compress_output bgzip \
                --fasta $FASTA \
                --fork ${task.cpus} \
                --vcf \
                --per_gene \
                --sift p \
        	--polyphen p \
	        --check_existing \
		--canonical \
		--buffer_site 5000	
	
	"""
}

workflow.onComplete {

  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

}

