#!/usr/bin/env nextflow

// MAIN PARAMETERS

if (!params.accessions && !params.fastq && !params.fasta){ exit 1, "'accessions', 'fasta' or 'fastq' parameter missing" }

if (params.fasta){
    if (params.fasta instanceof Boolean){exit 1, "'fasta' must be a path pattern. Provide value:'$params.fasta'"}

    IN_fasta_raw = Channel.fromPath(params.fasta).map{ it -> file(it).exists() ? [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it] : null }.ifEmpty { exit 1, "No fasta files provided with pattern:'${params.fasta}'" }

    IN_alignment_fasta = Channel.value(params.alignmentFasta)

    process RGI_FASTA{

        tag { sample_id }
        publishDir "results/rgi_fasta/"

        input:
        set sample_id, file(fasta) from  IN_fasta_raw
        val alignmetTool from IN_alignment_fasta

        output:
        file("*card_rgi*")

        script:
        """
        # Place card_rgi source in a read/write location for container
        mkdir card_temp && cp -r /opt/conda/lib/python3.6/site-packages/app/ card_temp
        export PYTHONPATH="\$(pwd)/card_temp/:\$PATH"
        rgi main --input_sequence ${fasta} --output_file ${sample_id}_card_rgi --input_type contig --alignment_tool ${alignmetTool} --low_quality --include_loose -d wgs --clean
        rgi parser -i ${sample_id}_card_rgi.json -o parsed_${sample_id}_card_rgi -t contig
        """

    }
    

} else {

    if (params.accessions){
        IN_accessions_raw = Channel.fromPath(params.accessions).ifEmpty { exit 1, "No accessions file provided with path:'${params.accessions}'" }

        process fasterqDump {

        tag { accession_id }
        maxRetries 1

        input:
        val accession_id from IN_accessions_raw.splitText(){ it.trim() }.filter{ it.trim() != "" }

        output:
        set val({ "$name" != "null" ? "$name" : "$accession_id" }), file("${accession_id}/*fq") optional true into IN_fastq_raw

        script:
        """
        {
            echo "Downloading the following accession: ${accession_id}"
            fasterq-dump ${accession_id} -e ${task.cpus} -p
            if [ ${params.compress_fastq} = true ]
            then
                echo "Compressing FastQ files..."
                if [ -f ${accession_id}_1.fastq ]
                then
                    pigz -p ${task.cpus} ${accession_id}_1.fastq ${accession_id}_2.fastq
                elif [ -f ${accession_id}_3.fastq ]
                then
                    echo "No paired end reads were found to compress."
                    pigz -p ${task.cpus} ${accession_id}_3.fastq
                else
                    echo "FastQ files weren't compressed. Check if FastQ files were downloaded."
                fi
            elsenex
                echo "FastQ files won't be compressed because compress_fastq options was set to: '${params.compress_fastq}.'"
            fi
        } || {
            # If exit code other than 0
            if [ \$? -eq 0 ]
            then
                echo "pass" > .status
            else
                echo "fail" > .status
                echo "Could not download accession $accession_id" > .fail
            fi
        }
        """

        }
    } else {

        if (!params.fastq){ exit 1, "'fastq' parameter missing"}
        // size: -1 -> allows for single and paired-end files to be passed through. Change if necessary
        IN_fastq_raw = Channel.fromFilePairs(params.fastq, size: -1).ifEmpty {
            exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

    }

    // PreProcessing
    process fastp {

        tag { sample_id }

        input:
        set sample_id, file(fastq_pair) from IN_fastq_raw

        output:
        set sample_id, file("*trim_*.fq.gz") into OUT_fastp

        script:
        """
        a=(${fastq_pair})
        if ((\${#a[@]} > 1));
        then
            fastp -i ${fastq_pair[0]} -o ${sample_id}_trim_1.fq.gz -I ${fastq_pair[1]} -O ${sample_id}_trim_2.fq.gz 
        else
            fastp -i ${fastq_pair} -o ${sample_id}_trim_1.fq.gz 
        fi
        """
    }

    IN_alignment_reads = Channel.value(params.alignmentReads)


    // RGI_BWT
    process RGI_BWT {

        tag { sample_id }
        publishDir "results/rgi_bwt/"
        
        input:
        set sample_id, file(fastq_pair) from OUT_fastp
        val alignmetTool from IN_alignment_reads

        output:
        "*_rgi_bwt*"

        script:
        """
        # Place card_rgi source in a read/write location for container
        mkdir card_temp && cp -r /opt/conda/lib/python3.6/site-packages/app/ card_temp
        export PYTHONPATH="\$(pwd)/card_temp/:\$PATH"
        rgi bwt --read_one ${fastq_pair[0]} --read_two ${fastq_pair[0]} --output_file ${sample_id}_rgi_bwt --aligner ${alignmetTool} --clean
        """
    }
}