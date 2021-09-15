![](images/ikmb_bfx_logo.png)

# IKMB VEP Pipeline

Run VEP on any number of VCFs

## Usage

To run this pipeline, you need nextflow and singularity. Then simply do:

`nextflow run ikmb/vep --assembly GRCh38 --vcfs /path/to/*.vcf.gz`

Currently, only the assembly version GRCh38 is supported.

