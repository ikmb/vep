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

// WORKFLOW starts here

process vep {

        label 'vep'

	publishDir "${OUTDIR}/DeepVariant/VEP", mode: 'copy'


	input:
        set file(vcf),file(vcf_index) from vcfs

	output:
        file(vcf_vep)
        file(vcf_alissa)

	script:
	vcf_vep = vcf.getBaseName() + ".vep.vcf"
	vcf_alissa = vcf.getBaseName() + ".vep2alissa.vcf"

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
                --fasta $FASTA \
                --fork 4 \
                --vcf \
                --per_gene \
                --sift p \
        	--polyphen p \
	        --check_existing \
		--canonical

		sed -i.bak 's/CADD_PHRED/CADD_phred/g' $vcf_vep

	"""
}


// Header log info
log.info "========================================="
log.info "VEP pipeline v${params.version}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "-----------------------------------------"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${run_name}"
if (workflow.containerEngine) {
	log.info "Container engine:		${workflow.containerEngine}"
}
log.info "========================================="


workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['Samples'] = params.samples
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['kit'] = TARGETS
  email_fields['assembly'] = FASTA
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "Diagnostic exome analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[IKMB ExoSeq] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[IKMB ExoSeq] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}

