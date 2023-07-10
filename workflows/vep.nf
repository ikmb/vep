fasta = params.genomes[params.assembly].fasta
fai = params.genomes[params.assembly].fai

ch_fasta = Channel.from([fasta,fai])

include { ENSEMBL_VEP } from '../modules/ensembl/vep.nf'
include { VEP2XLS } from '../modules/vep2xls'
include { MULTIQC } from '../modules/multiqc'
include { SOFTWARE_VERSIONS } from '../modules/software_versions'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions/main'

Channel.fromPath(params.vcf).map { f ->
    [[
        id: f.getSimpleName()
    ],
    f,file(f + ".tbi", checkIfExists: true)
    ]
}.set { ch_vcfs }

ch_versions = Channel.from([])
multiqc_files = Channel.from([])

log.info "Running with options: ${params.vep_options}"

workflow VEP {

    main:

    ENSEMBL_VEP(
        ch_vcfs,
        ch_fasta.collect()
    )

    ch_versions = ch_versions.mix(ENSEMBL_VEP.out.versions)

    VEP2XLS(
        ENSEMBL_VEP.out.vcf
    )

}
