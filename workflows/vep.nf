fasta = params.genomes[params.assembly].fasta
fai = params.genomes[params.assembly].fai

ch_fasta = Channel.from([fasta,fai])

include { ENSEMBL_VEP } from '../modules/ensembl/vep.nf'
include { MULTIQC } from '../modules/multiqc'
include { SOFTWARE_VERSIONS } from '../modules/software_versions'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

Channel.fromPath(params.vcf).map { f ->
    [[
        id: f.getSimpleName()
    ],
    f,file(f + ".tbi", checkIfExists: true)
    ]
}.set { ch_vcfs }

ch_versions = Channel.from([])
multiqc_files = Channel.from([])

workflow VEP {

    take:
    vcf

    main:

    ENSEMBL_VEP(
        vcf,
        ch_fasta.collect()
    )

    ch_versions = ch_versions.mix(ENSEMBL_VEP.out.versions)
    multiqc_files = multiqc_files.mix(FASTP.out.json)

    SOFTWARE_VERSIONS(
        ch_versions.collect()
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    MULTIQC(
        multiqc_files
    )

    emit:
    qc = MULTIQC.out.report

}
