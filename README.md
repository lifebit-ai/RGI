# RGI
Antimicrobial resistance gene identification through CARD's Resistance Gene Identifier

## Rationale
This workflow performs *in silico* screening of **antimicrobial resistance genes** in genomic data through [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi). This relies on the [Comprehensive Antibiotic Resistance (CARD) Database](https://card.mcmaster.ca/). 

Two RGI modes are implemented in this workflow: the `main` and the `bwt`. The first recieves as input already pre-assembled FASTA sequences, where the second aligns sequence reads (FASTQ)to CARD to make the prediction. The main algorithms implemented in each are:
- main: BLAST or DIAMOND
- bwt: Bowtie2 or BWA

## Usage

