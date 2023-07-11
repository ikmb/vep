process VEP2XLS {

    label 'medium_serial'

    container 'ikmb/vep:1.0'

    publishDir "${params.outdir}/VEP", mode: 'copy'

    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path(sheet), emit: xlsx

    script:
    sheet = vcf.getBaseName() + ".xlsx"

    """
    zcat $vcf | vep2xls_fast.rb  -o $sheet
    """

}

