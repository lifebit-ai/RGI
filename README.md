# RGI
Antimicrobial resistance gene identification through CARD's Resistance Gene Identifier

## Rationale
This workflow performs *in silico* screening of **antimicrobial resistance genes** in genomic data through [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi). This relies on the [Comprehensive Antibiotic Resistance (CARD) Database](https://card.mcmaster.ca/). 

## Implementation
Two RGI modes are implemented in this workflow: the `fasta` and the `fastq`. The first recieves as input already pre-assembled FASTA sequences, where the second aligns sequence reads (FASTQ) to CARD to make the prediction. The main algorithms implemented in each are:
- fasta: BLAST or DIAMOND (default: DIAMOND)
- fastq: Bowtie2, bwa or kma (default: kma)

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
    --accessions        Path file with accessions, one perline.
                        Input type: path (default: )

    FastqDump Options:
    --compress_fastq    If true, downloaded fastqs will be compressed.
                        Input tupe: boolean (default: true)

    RGI Options:
    --alignmentFasta    Specify alignment tool. Options: {DIAMOND,BLAST}
                        Input type: string (default: DIAMOND)
    --alignmentReads    Specify alignment tool. Options: {bowtie2,bwa,kma}
                        Input type: numeric (default: kma)

## Basic run command example
    nextflow main.nf --fastq="test_data/SRR8886136_{1,2}*"