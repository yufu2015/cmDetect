#ctDNA Mutation Detection Workflow (cmDetect)
cmDetect is a method for the systematic identification of ctDNA mutations in a large set of whole-exome sequencing of tumor and corresponding peripheral whole-blood samples
Contact: yu.fu@gustaveroussy.fr

## Installation
no need for installation

## External Prerequisites:

cmDetect requires the following software:
GATK (version 3.4.0 or later)
https://www.broadinstitute.org/gatk/download/
SnpEff (version 4.1 or later)
http://snpeff.sourceforge.net/download.html
Samtools (version 0.1.19 or later)
http://sourceforge.net/projects/samtools/files/samtools
R (version 3.2.1) 
https://cran.r-project.org/

Note: the version of the tools provided here are the ones used in our analysis, earlier version might work as well except for samtools


## Usage

A detailed cmDetect workflow on a single procesor is described in pipeline.sh in the script folder.
If you wish to run certain steps (variant calling using GATK, annotation etc.) parallely, please adapt all the scripts accordingly (set environment variables etc.).
It is recommended to run the workflow step by step to ensure that each one is properly finished for all the samples.

## Help
A log file with all the error messages can be found under the cmDetect directory
 


