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
--vcfs			       List of VCFs to process
--email 		       Email address to send reports to (enclosed in '')
--sites			       Remove individual genotype data to speed up processing of very large files
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
	log.info "!!!Please consider using a more recent genome version!!!"
}

if (!params.vep_cache_dir || !params.vep_plugin_dir) {
	exit 1, "Missing VEP cache and/or plugin directory..."
}

dbNSFP_DB = file(params.genomes[ params.assembly ].dbnsfp_db)
if (!dbNSFP_DB.exists()) {
	exit 1, "Could not find the specified dbNSFP database..."
}

dbscSNV_DB = file(params.genomes[ params.assembly].dbscsnv_db)
if (!dbscSNV_DB.exists()) {
	dbscSNV_DB = file(params.dbscsnv_db)
	exit 1, "Could not find the specified dbscSNV database..."
}

CADD_SNPS = file(params.genomes[ params.assembly].cadd_snps)
CADD_INDELS = file(params.genomes[ params.assembly].cadd_indels)
if ( !CADD_SNPS.exists() || !CADD_INDELS.exists() ) {
	exit 1, "Missing CADD SNPs and/or Indel references..."
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

if (params.sites) {

	process vcf_sites_only {

		label 'gatk'

		publishDir "${params.outdir}/VCFS_SITES_ONLY", mode: 'copy'

		input:
		set file(vcf),file(vcf_index) from vcfs

		output:
		set file(vcf_sites),file(vcf_sites_index) into vcf_sites

		script:
	
		vcf_sites = vcf.getBaseName() + ".sites.vcf.gz"
		vcf_sites_index = vcf_sites + ".tbi"

		"""
			gatk SelectVariants -V $vcf --sites-only-vcf-output -O $vcf_sites -OVI
		"""
	}

} else {

	vcf_sites = vcfs

}

process vep {

        label 'vep'

	publishDir "${params.outdir}/VEP", mode: 'copy'

	input:
        set file(vcf),file(vcf_index) from vcf_sites

	output:
        file(vcf_vep)

	script:
	vcf_vep = vcf.getBaseName() + ".vep.txt.gz"
	options = ""
	if (params.revel) {
		options = "--plugin REVEL,${params.revel}"
	}
	if (params.refseq) {
		options = options + " --refseq"
	}
	"""
		export PERL5LIB=${params.vep_plugin_dir}

                vep --offline \
                --cache \
                --dir ${params.vep_cache_dir} \
                --species homo_sapiens \
                --assembly $params.assembly \
                -i $vcf \
                --format vcf \
		--tab \
		--hgvs \
                -o $vcf_vep --dir_plugins ${params.vep_plugin_dir} \
                --plugin dbNSFP,$dbNSFP_DB,${params.dbnsfp_fields} \
                --plugin dbscSNV,$dbscSNV_DB \
                --plugin CADD,${params.cadd_snps},${params.cadd_indels} \
                --plugin ExACpLI \
		--plugin LoFtool \
		--plugin UTRannotator \
		--compress_output bgzip \
                --fasta $FASTA \
                --fork ${task.cpus} \
                --per_gene \
                --sift p \
        	--polyphen p \
	        --check_existing \
		--canonical \
		$options
	
	"""
}

workflow.onComplete {

  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

}

