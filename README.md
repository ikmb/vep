![](images/ikmb_bfx_logo.png)

# IKMB VEP Pipeline

Run VEP on any number of VCF files. This pipeline only works on the IKMB MEDCluster due to the inclusion of VEP plugins and flatfile databases. 

## Usage

To run this pipeline, you need nextflow and singularity. Then simply do:

`nextflow run ikmb/vep --assembly GRCh38 --vcfs '/path/to/*.vcf.gz'`

Currently, only the assembly version GRCh38 is supported.

## Options

### `--sites` [ true | false (default) ]

Whether to remove individual genotype columns from the VCF first (useful for very large VCF files)

### `--revel`[ true | false (default) ]

Whether to produce REVEL scores as well
