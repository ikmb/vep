process ENSEMBL_VEP {

    publishDir "${params.outdir}/VEP", mode: 'copy'
    
    label 'medium_parallel'

    container 'quay.io/biocontainers/ensembl-vep:109.3--pl5321h2a3209d_1'

    input:
    tuple val(meta),path(v),path(t)
    tuple path(fasta),path(fai)

    output:
    tuple val(meta),path(vep_vcf), emit: vcf
    path("versions.yml"), emit: versions

    script:

    """
    export PERL5LIB=${params.vep_plugin_dir}

    vep --offline \
        --cache \
        --dir ${params.vep_cache_dir} \
        --species homo_sapiens \
        --assembly $params.assembly \
        -i $vcf \
        --format vcf \
        --hgvs \
        -o $vcf_vep --dir_plugins ${params.vep_plugin_dir} \
        --plugin dbNSFP,${params.dbnsfp_db},${params.dbnsfp_fields} \
        --plugin dbscSNV,${params.dbscsnv_db} \
        --plugin CADD,${params.cadd_snps},${params.cadd_indels} \
        --plugin ExACpLI \
        --plugin UTRannotator \
        --plugin Mastermind,${params.vep_mastermind} \
        --plugin SpliceAI,${params.spliceai_fields} \
        --af_gnomadg \
        --compress_output bgzip \
        --fasta $fasta \
        --fork ${task.cpus} \
        --per_gene \
        --sift p \
        --polyphen p \
        --check_existing \
        --canonical \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep \$( vep -help 2>&1 | grep "ensembl-vep" | sed -e "s/.*: //g" )
    END_VERSIONS

    """	

}


