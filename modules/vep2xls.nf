process VEP2XLS {

    label 'medium_serial'

    container 'ikmb/exome-seq:5.2'

    publishDir "${params.outdir}/VEP", mode: 'copy'

    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path(sheet), emit: xlsx

    script:
    sheet = vcf.getBaseName() + ".xlsx"

    """
    zcat $vcf > variants.vcf
    vep2xls.rb -i variants.vcf -o $sheet
    rm variants.vcf
    """

}

