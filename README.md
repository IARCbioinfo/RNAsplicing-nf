# RNAsplicing-nf
## Nextflow pipeline to perform RNA splicing analyses

[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/rnasplicing-nf/)
[![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)

![Workflow representation](template-nf.png)

## Description
Runs the SUPPA2 software after recommended data processing of fastq files.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- trimgalore
- salmon
- SUPPA2

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input
  | Type      | Description     |
  |-----------|---------------|
  | input_folder    | Folder containing fastq files to be aligned |
  | input_file    | Input file (tab-separated values) with 4 columns: ID (sample name), RG (read group), pair1 (first fastq pair file), and pair2 (second fastq pair file) |

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
 | index       |  salmon_index  |  Reference fasta file (with index) for splice junction trimming and base recalibration |
 | gtf          |  GRCh38_gencode_annotation.gtf  |   Path to annotation gtf file |


  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| output_folder    | RNAsplicing-nf_output |      Output folder'
| cpu                  | 4 | Number of cpu used by salmon (default: 4).'
| mem                             | 10| Size of memory used for mapping (in GB) '
| fastq_ext                       | fq.gz | Extension of fastq files '
| suffix1                         | 1| Suffix of fastq files 1'
| suffix2                         | 2|                  Suffix of fastq files 2'
| mem_QC                         | 2|                  Size of memory used for QC and cutadapt (in GB)'
| cpu_trim                       |  15|                  Number of cpu used by cutadapt '
| index                           | salmon_index |                  Path to salmon index folder'
| suppa_folder                    | SUPPA |                  Path to SUPPA folder (cloned from github) containing suppa.py and folders scripts.'
| ioe    | null  |                    Path to SUPPA ioe file. If not specified, will generate the file.'


## Usage
  ```
 nextflow run iarcbioinfo/RNAsplicing-nf --input_file input.tsv --index salmon_index --gtf annot.gtf 
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | QC/adapter_trimming    | Fastqc reports from trim_galore |
  | quantification   | Transcript expression quantification from salmon |
  | results/*ioe and ioi files | SUPPA2 ioe and ioi files used for quantification of psi by events and isoforms  |
  | results/psiPerEvent   | Files with psi per event for each sample |
  | results/psiPerIsoform | Files with psi per isoform for each sample |


## Detailed description
Nextflow pipeline that implements the [SUPPA2 tutorial](https://github.com/comprna/SUPPA/wiki/SUPPA2-tutorial), first running salmon for quantification, then postprocessing the salmon outputs, creating the SUPPA2 ioe and ioi files if necessary, and computing psi per event and isoform for each sample.

## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/RNAsplicing-nf/blob/master/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Nicolas Alcala*    |            alcalan@iarc.who.int | Developer to contact for support |
  | Ricardo Blazquez Encinas-Rey    |            | Developer |

## References 
* Trincado JL, Entizne JC, Hysenaj G, Singh B, Skalic M, Elliott DJ, Eyras E. [SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions](https://www.ncbi.nlm.nih.gov/pubmed/29571299). Genome Biol. 2018 Mar 23;19(1):40.
* Alamancos GP, Pagès A, Trincado JL, Bellora N, Eyras E. [Leveraging transcript quantification for fast computation of alternative splicing profiles](https://www.ncbi.nlm.nih.gov/pubmed/26179515). RNA. 2015 Sep;21(9):1521-31.
