# Usage and Parameters

## Requirements

This workflow requires at least `4` cpus and `8GB` of memory.

## Usage

The typical command for running the pipeline is as follows:

    nextflow run main.nf [Options]

    Inputs Options:
    --fastq             Path expression to fastq files.
                        Input type: path (default: )
    --fasta             Path expression to fasta files.
                        Input type: string (default: )
    --accessions        Path file with accessions, one per line.
                        Input type: path (default: )

    FastqDump Options:
    --compress_fastq    If true, downloaded fastqs will be compressed.
                        Input tuple: boolean (default: true)

    RGI Options:
    --alignmentFasta    Specify alignment tool. Options: {DIAMOND,BLAST}
                        Input type: string (default: DIAMOND)
    --alignmentReads    Specify alignment tool. Options: {bowtie2,bwa,kma}
                        Input type: numeric (default: kma)

## Basic run command example

### fasta

    nextflow run main.nf --fasta="*.fasta"

### fastq

    nextflow run main.nf --fastq="test_data/SRR8886136_{1,2}*"

## Test run

    nextflow run main.nf --profile test
