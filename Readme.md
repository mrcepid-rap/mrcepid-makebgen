# MakeBGEN (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Background](#background)
  * [Dependencies](#dependencies)
    + [Docker](#docker)
    + [Resource Files](#resource-files)
- [Methodology](#methodology)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs)
  * [Outputs](#outputs)
  * [Command line example](#command-line-example)

## Introduction

This applet converts filtered and annotated individual VCFs on a single chromosome into two files:

1. A [.bgen format](https://www.well.ox.ac.uk/~gav/bgen_format/) file containing all genotypes a matrix of samples x variants for every participant in the UKBB
2. A tab-delimited .tsv.gz file with annotations for all variants in the `.bgen` format file.

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```commandline
dx describe file-1234567890ABCDEFGHIJKLMN
```

**Note:** This README pertains to data included as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk)

### Background

VCF and BCF format files take up a comparitively large amount of space due to a large number of FORMAT fields attached to
each genotype. Therefore, we have decided to use .bgen format to store genotype information after quality control and 
annotation due to it's relatively small size, accompanied by a tab-delimited .tsv file storing the information per-variant
normally provided as INFO fields in the VCF/BCF spec. This allows for both rapid conversion of .bgen files into other
formats for downstream burden testing, and relative quick I/O when processing variant annotations when selecting variants
for burden testing.

Please see the [.bgen format specification](https://www.well.ox.ac.uk/~gav/bgen_format/) for more information on these files.
For the purposes of compatibility with several other pipelines, we use the following .bgen format specifications:
* version 1.2
* ref-last – plink2 default (**Note:** UK Biobank imputation bgen files are 'ref-first', which may lead to confusion when adapting code from this repository and others!)
* 8 bit precision – for compatibility with most burden testing tools

For information on the fields included in the .tsv annotation, please see our documentation for [mrcepid-filterbcf](https://github.com/mrcepid-rap/mrcepid-filterbcf).

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/burdentesting.Dockerfile

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually. For more details on how to build a Docker image for use on the UKBiobank RAP, please see:

https://github.com/mrcepid-rap#docker-images

In brief, the primary **bioinformatics software** dependencies required by this Applet (and provided in the associated Docker image)
are:

* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/index.html)
* [bgenix](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)
* [plink2](https://www.cog-genomics.org/plink/2.0/)

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

#### Resource Files

In order for the applet to safely decide how to create BGEN chunks (with chunk-ends being in non-exonic regions), we need to create a dictionary of the genome in 
terms of gene coordinates. You can create this dictionary by running the following command from the root of the repository:

```bash
python scripts/gene_dict.py
```

This will create a file called `final_dict_public.json` in the home directory of the repository. You should upload this file to DNA Nexus, and use the file-dxid as part 
of the applet for the `igene_dict` parameter.

## Methodology

This applet is step 3 (mrc-makebgen) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
RAP at the MRC Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.v3.png)

The app has two primary steps:

1. Convert all provided filtered and annotated BCFs for one chromosome into bgen1.2 format with `plink2`:

```shell
plink2 --threads 2 --memory 10000 \
                  --bcf chr1_chunk1.filtered_annotated.bcf  \
                  --export bgen-1.2 'bits='8 'sample-v2'  \
                  --vcf-half-call r  \
                  --out filtered_annotated \
                  --set-all-var-ids '@:#:\$r:\$a' \
                  --new-id-max-allele-len 100
```

2. Merge the resulting `.bgen` files into a single whole-chromosome `.bgen` with the `cat-bgen` tool provided with the BGEN library:

```shell
cat-bgen -g chr1_chunk1.filtered_annotated.bgen -g chr1_chunk1.filtered_annotated.bgen -g chr1_chunk1.filtered_annotated.bgen \
        -og chr1.filtered_annotated.bgen
```

As part of step (2), the corresponding VEP annotations are merged as a DataFrame from the `pandas` package for python3.

## Running on DNANexus

### Inputs

| input           | description                                                                                                                           |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------|
| output_prefix   | the prefix for output files                                                                                                           |
| coordinate_file | a gzipped .tsv file of file coordinates and DNANexus file-ids to merge. See below for the format of this file.                        |
| gene_dict       | a file containing a dictionary of gene coordinates for the genome. This is used to determine chunk ends. **[final_dict_public.json]** |
| size_of_bgen    | desired size of the final merged BGEN file in megabases (e.g. 10 would indicate a final 10Mb BGEN file). 3Mb is the minimum.          |
| make_bcf        | Make a bcf file with identical sites/genotypes to the output bgen? **[False]**                                                        |

Each job from [mrcepid-filterbcf](https://github.com/mrcepid-rap/mrcepid-filterbcf) will produce a single coordinates file
for the files processed in that job. The user must merge these files on their own to provide the required input for this 
applet. The easiest way to do this is to launch a [Cloud Workstation](https://documentation.dnanexus.com/developer/cloud-workstations/cloud-workstation)
on DNANexus, download all coordinates files. 

**MASSIVE NOTE**: It is 1000% essential that the coordinate file is sorted as demonstrated below. Otherwise duplicate
variants spread across multiple bcfs will be included in the final bgen!

```shell
# Make a directory to do processing and download files
mkdir coords & cd coords
dx download -a --no-progress --lightweight filtered_annotated_vcfs/coordinates_*.tsv

## Do some quick QC to make sure we have 4722 UNIQUE files:
# Should give 4722 as output
cat *.tsv | awk '{print $4}' | sort | uniq | wc -l 

# Make a header file...
echo -e '#chrom\tstart\tstop\tfileprefix\tbcf_dxpy\tvep_dxpy' > coord_header.tsv
# Concatenate everything together:
cat coord_header.tsv > bcf_coordinates.tsv; cat *.tsv | sort -k 1,1 -k 2,2n coordinates_*.tsv >> bcf_coordinates.tsv
bgzip 

# Upload
dx upload --destination filtered_annotated_vcfs/ bcf_coordinates.tsv.gz
```

This file should have the following columns **WITH** a header:

| column name | description                                                                                                                                                                                                                         |
|-------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chrom       | chromosome, in `chr1`-like format                                                                                                                                                                                                   |
| start       | position of the first variant in this bcf                                                                                                                                                                                           |
| end         | position of the last variant in this bcf                                                                                                                                                                                            |
| vcf_prefix  | prefix string for this bcf. Will be something that matches the regex to `ukb\d+_c\d_b\d*_v1_chunk\d`, where the numbers (`\d`) in order are: UKBB field for current WES release, chromosome, original vcf slice, current vcf slice. |
| output_bcf  | DNANexus file ID (like `file-1234567890ABCDEFGHIJKLMN`) for the `.bcf` indicated by vcf_prefix                                                                                                                                      |
| output_vep  | DNANexus file ID (like `file-1234567890ABCDEFGHIJKLMN`) for the corresponding annotations .tsv for the `.bcf` indicated by vcf_prefix.                                                                                              |

### Outputs

This applet generates two primary outputs per chromosome:

1. A [.bgen format](https://www.well.ox.ac.uk/~gav/bgen_format/) file
2. A [bgzipped](http://www.htslib.org/doc/bgzip.html) tab-delimited .tsv.gz file with annotations for all variants in the `.bgen` format file.

These files are accompanied by various indices as indicated below:

| output   | description                                                                                          |
|----------|------------------------------------------------------------------------------------------------------|
| bgen     | the `.bgen` format file                                                                              |
| index    | The corresponding .bgi index for the .bgen format file                                               |
| sample   | A sample file for the .bgen                                                                          |
| vep      | A bgzipped `.tsv.gz` of annotations for all variants included in the corresponding `.bgen` file      |
| vep_idx  | A [tabix index](http://www.htslib.org/doc/tabix.html) for the `.tsv.gz` allowing query by coordinate |
| bcf      | _Optional_: A concatenated BCF file with identical variants/genotypes to bgen                        |
| bcf_indx | _Optional_: The index for the concatenated BCF file                                                  |



### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our 
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

Running this command is fairly straightforward using the DNANexus SDK toolkit:

```shell
dx run mrcepid-makebgen --priority low --destination filtered_bgen/ -ichromosome='chr1' -icoordinate_file=file-1234567890ABCDEFGHIJKLMN
```

Some notes here regarding execution:

1. For ease of execution, I prefer using the file hash. This is mostly because DNANexus has put lots of spaces in their 
   filepaths, AND it is easier to programmatically access many files at once using hashes.

2. Outputs are automatically named based on the prefix of the input vcf full path (this is regardless of if you use hash or full path). So 
   the primary bgen output for the above command-line will be `chr1.filtered.bgen`. All outputs 
   will be named using a similar convention.

3. I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for building the app (at `dxapp.json`) 
   so setting an instance type is unnecessary. This current default is for a mem3_ssd1_v2_x16 instance (16 CPUs, 128 Gb RAM, 600Gb storage). This
   instance prioritises more RAM over other types of instances. **Please note** that this applet is set up for the 
   parallelisation of many files. To run one file, one needs much less memory. If necessary to adjust compute resources,
   one can provide a flag like `--instance-type mem3_ssd1_v2_x8` to `dx run`.