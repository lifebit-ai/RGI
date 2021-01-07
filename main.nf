#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fasta sample.fasta [Options]
    
    Inputs Options:
    --fastq             Path expression to fastq files.
                        Input type: path (default: $params.fastq)
    --fasta             Path expression to fasta files.
                        Input type: string (default: $params.fasta)
    --accessions        Path file with accessions, one perline.
                        Input type: path (default: $params.accessions)
    
    FastqDump Options:
    --compress_fastq    If true, downloaded fastqs will be compressed. 
                        Input tupe: boolean (default: $params.compress_fastq)
    
    RGI Options: 
    --alignmentFasta    Specify alignment tool. Options: {DIAMOND,BLAST}
                        Input type: string (default: $params.alignmentFasta)
    --alignmentReads    Specify alignment tool. Options: {bowtie2,bwa,kma}
                        Input type: numeric (default: $params.alignmentReads)
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// MAIN PARAMETERS

if (!params.accessions && !params.fastq && !params.fasta){ exit 1, "'accessions', 'fasta' or 'fastq' parameter missing" }

if (params.fasta){
    if (params.fasta instanceof Boolean){exit 1, "'fasta' must be a path pattern. Provide value:'$params.fasta'"}

    IN_fasta_raw = Channel.fromPath(params.fasta).map{ it -> file(it).exists() ? [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it] : null }.ifEmpty { exit 1, "No fasta files provided with pattern:'${params.fasta}'" }

    IN_alignment_fasta = Channel.value(params.alignmentFasta)

    process RGI_FASTA{

        tag { sample_id }
        publishDir "results/rgi_fasta/", pattern: "*.txt"

        input:
        set sample_id, file(fasta) from  IN_fasta_raw
        val alignmetTool from IN_alignment_fasta

        output:
        file("*card_rgi.json") into CARD_HEATMAP
        file("*_card_rgi_parsed-count-hits.json") into OUT_RGI_FASTA

        script:
        """
        # Place card_rgi source in a read/write location for container
        mkdir card_temp && cp -r /opt/conda/lib/python3.6/site-packages/app/ card_temp
        export PYTHONPATH="\$(pwd)/card_temp/:\$PATH"
        
        rgi main --input_sequence ${fasta} --output_file ${sample_id}_card_rgi --input_type contig --alignment_tool ${alignmetTool} --low_quality --include_loose -d wgs --clean -n $task.cpus
        
        rgi parser -i ${sample_id}_card_rgi.json -o ${sample_id}_card_rgi_parsed --include_loose -t contig
        
        #clean up work dir, if it exists
        [[ -d card_temp ]] && rm -r card_temp
        """
    }

    process PROCESS_RGI_HEATMAP {

        publishDir "results/MultiQC/", pattern: "*.png", mode: "copy"

        input:
        file(JSON_HITS) from CARD_HEATMAP.collect()

        output:
        file("results_hits.csv") into DF_TABLE_HITS
        file("card_hits.csv") into DF_HEATMAP

        script:
        template "process_json_hits.py"
    }

    process PROCESS_RGI_FASTA {

        input:
        file(JSON_FILES) from OUT_RGI_FASTA.collect()

        output:
        file("results_summary.csv") into OUT_STATS_SUMMARY

        script:
        template "parse_rgi_json.py"
    }


} else {

    if (params.accessions){
        IN_accessions_raw = Channel.fromPath(params.accessions).ifEmpty { exit 1, "No accessions file provided with path:'${params.accessions}'" }

        process fasterqDump {

        tag { accession_id }
        publishDir "reads/${accession_id}/", pattern: "*.fastq*"
        maxRetries 1

        input:
        val accession_id from IN_accessions_raw.splitText(){ it.trim() }.filter{ it.trim() != "" }
        
        output:
        set accession_id, file("*.fastq*") optional true into IN_fastq_raw

        script:
        """
        echo "Downloading the following accession: ${accession_id}"
        fasterq-dump ${accession_id} -e ${task.cpus}
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
        else
            echo "FastQ files won't be compressed because compress_fastq options was set to: '${params.compress_fastq}.'"
        fi
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
        file("*_mapping_data.txt")
        file("*mapping_stats.txt")
        file("*.json") into OUT_RGI_JSON_BWT

        script:
        """
        # Place card_rgi source in a read/write location for container
        mkdir card_temp && cp -r /opt/conda/lib/python3.6/site-packages/app/ card_temp
        export PYTHONPATH="\$(pwd)/card_temp/:\$PATH"
        
        rgi bwt --read_one ${fastq_pair[0]} --read_two ${fastq_pair[0]} --output_file ${sample_id}_rgi_bwt --aligner ${alignmetTool} --clean
        """
    }

    process PROCESS_RGI_BWT {


        input:
        file JSON_FILES from OUT_RGI_JSON_BWT.collect()

        output:
        file("results_summary.csv") into OUT_STATS_SUMMARY
        file("results_hits.csv") into DF_TABLE_HITS
        file("card_hits.csv") into DF_HEATMAP

        script:
        template "process_rgi_bwt.py"

    }
}

/*
process REPORT {

    publishDir "results/MultiQC/", mode: "copy"

    input:
    file summary_html from OUT_STATS_SUMMARY
    file html_hit_table from HTML_TABLE_HITS
    file heatmap from PNG_HEATMAP
    
    output:
    file "*report.html" into OUT_REPORT

    script:
    template "generate_report.py"
}
*/

process REPORT {

    publishDir "results/MultiQC/", mode: "copy"

    input:
    file summary_df from OUT_STATS_SUMMARY
    file hit_table from DF_TABLE_HITS
    file heatmap_df from DF_HEATMAP
    
    output:
    file "*report.html" into OUT_REPORT

    script:
    template "report.py"
}

workflow.onComplete {

  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}