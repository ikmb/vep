#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
VEP Pipeline
===============================

This Pipeline performs variant effect prediction using VEP

### Homepage / git
git@github.com:ikmb/pipeline.git

**/

// Pipeline version

params.version = workflow.manifest.version

def summary = [:]

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

WorkflowMain.initialise(workflow, params, log)
WorkflowVep.initialise( params, log)

include { VEP } from './workflows/vep'

workflow {

    VEP()

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

}

