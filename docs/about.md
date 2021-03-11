# RGI

Antimicrobial resistance gene identification through CARD's Resistance Gene Identifier


## Rationale

This workflow performs *in silico* screening of **antimicrobial resistance genes** in genomic data through [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi). This relies on the [Comprehensive Antibiotic Resistance (CARD) Database](https://card.mcmaster.ca/). 


## Implementation
Two RGI modes are implemented in this workflow: the `fasta` and the `fastq`. 

### FASTA

This mode input receives already pre-assembled sequences in FASTA format. Two alignment tools are available for the detection of AMR:

- [DIAMOND](https://github.com/bbuchfink/diamond) (**default**): a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi): local alignment tool to find regions of similarity between biological sequences.

### FASTQ

This mode aligns sequence reads (FASTQ) to CARD Database to make the prediction. The main alignment tools  implemented are:

- [kma](https://bitbucket.org/genomicepidemiology/kma/src/master/) (**default**): a mapper of raw reads directly against redundant databases, in an ultra-fast manner using k-mer seed and extension.
- [bwa](http://bio-bwa.sourceforge.net/): Implements BWA-MEM, generally recommended for high-quality short read data for higher accuracy and better performance 
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): aligns sequencing reads to long reference sequences, particularly good at with reads of about 50 up to 100 basepairs. Implements the very-sensitive-local algorithm.

