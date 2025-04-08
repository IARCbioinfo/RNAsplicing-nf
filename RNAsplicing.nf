#! /usr/bin/env nextflow

// Copyright (C) 2025 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help         = null

params.input_folder = null
params.input_file   = null

params.cpu          = 4
params.cpu_trim     = 15
params.mem          = 50
params.mem_QC       = 2
params.fastq_ext    = "fq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.output_folder= "."
params.index          = "ref.ind"


log.info ""
log.info "--------------------------------------------------------"
log.info "  RNAsplicing-nf v1.0: RNA splicing analysis         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/RNAsplicing-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '--input_folder   FOLDER                 Folder containing fastq files to be aligned.'
    log.info '--input_file     STRING                 Input file (tab-separated values) with 4 columns: SM (sample name), RG (read group), pair1 (first fastq pair file), and pair2 (second fastq pair file)'
    log.info '    --index            FILE                   Reference fasta file (with index) for splice junction trimming and base recalibration.'
    log.info ""
    log.info "Optional arguments:"
    log.info '--output_folder                   STRING                 Output folder (default: .).'
    log.info '--cpu                             INTEGER                Number of cpu used by bwa mem and sambamba (default: 4).'
    log.info '--mem                             INTEGER                Size of memory used for mapping (in GB) (default: 50).'
    log.info '--fastq_ext                       STRING                 Extension of fastq files (default : fq.gz)'
    log.info '--suffix1                         STRING                 Suffix of fastq files 1 (default : _1)'
    log.info '--suffix2                         STRING                 Suffix of fastq files 2 (default : _2)'
    log.info '--mem_QC                         INTEGER                    Size of memory used for QC and cutadapt (in GB) (default: 2).'
    log.info '--cpu_trim                       INTEGER                    Number of cpu used by cutadapt (default: 15).'
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
log.info "input_folder   = ${params.input_folder}"
log.info "input_file     = ${params.input_file}"
log.info "index          = ${params.index}"
log.info "fastq_ext      = ${params.fastq_ext}"
log.info "suffix1        = ${params.suffix1}"
log.info "suffix2        = ${params.suffix2}"
log.info "output_folder  = ${params.output_folder}"
log.info "cpu            = ${params.cpu}"
log.info "mem_QC         = ${params.mem_QC}"
log.info "cpu_trim       = ${params.cpu_trim}"
}


process trim{
    cpus params.cpu_trim
    memory params.mem_QC+'GB'
    tag { file_tag +rg }
	    
    input:
        tuple val(file_tag), val(rg), path(pair1), path(pair2)
	    
    output:
        tuple val(file_tag), val(rg) , path("${file_tag}${rg}*val_1.fq.gz"), path("${file_tag}${rg}*val_2.fq.gz")  
	    path("*_fastqc.zip")
	    path("*trimming_report.txt")
        
    publishDir "${params.output_folder}/QC/adapter_trimming", mode: 'copy', pattern: '{*report.txt,*fastqc.zip}'
	    
    shell:
	    cpu_tg = params.cpu_trim -1
	    cpu_tg2 = cpu_tg.div(3.5)
	    cpu_tg3 = Math.round(Math.ceil(cpu_tg2))
        if(suffix2){
            pairs="${pair1} ${pair2}"
            opts="--paired "
        }else{
            pairs="${pair1}"
            opts=" "
        }
        '''
	    trim_galore !{opts} --fastqc --gzip --basename !{file_tag}!{rg} -j !{cpu_tg3} !{pairs}
        if [ -L NO_fastq2 ]
            mv !{file_tag}!{rg}_trimmed.fq.gz !{file_tag}!{rg}_val_1.fq.gz
            then touch !{file_tag}!{rg}_val_2.fq.gz 
        fi 
        '''
}

process quant{
    cpus params.cpu
    memory params.mem+'GB'
    tag {ID}

    input:
        path index
        tuple val(file_tag), val(rg) , path(file1), path(file2)  
	    
    output:
        tuple val(file_tag), path "file_tag"
    publishDir "${params.output_folder}/quantification", mode: "copy"
    
    script:
    """
    salmon quant -i ${index} -l A --gcBias --validateMappings -1 ${file1} -2 ${file2} -p 2 -o ${file_tag}
    """
}


process SUPPA2{
    cpus params.cpu
    memory params.mem+'GB'
    tag {file_tag}

    input:
        tuple val(file_tag), path salmon_folder

    output:
        path "*"
    publishDir "${params.output_folder}/results", mode: "copy"
    
    script:
    """
    python $projectDir/bin/suppa2_run.py $salmon_folder
    """
}

workflow{
    if(params.input_file){
	bams = Channel.fromPath("${params.input_file}")
     	          .splitCsv( header: true, sep: '\t', strip: true )
	       .map { row -> [ row.ID , row.RG, file(row.input_folder+"*1.fastq.gz") , file(row.matrix_folder+"*2.fastq.gz") ] }
	       .view()
    }else{
        fastqs     = Channel.fromFilePairs( params.fastq+'{*1.fastq.gz,*2.fastq.gz}' ).view()
    }
    index     = file(params.index)
    
    output_trim_ch = trim(fastqs)
    output_quant_ch = quant(output_trim_ch,index)
    SUPPA2(output_quant_ch[0])
}