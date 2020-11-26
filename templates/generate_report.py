#!/usr/bin/env python3

JSON_REPORTS = "$json_files".split()

html_template = """<html>
<head></head>
<body>
    
    <h1> RGI Workflow </h1>
    <h3>About</h3>

    <p> 
    A pipeline to perform antimicrobial resistance detection with [POLYSOLVER](https://software.broadinstitute.org/cancer/cga/polysolver).

Input can be of three type -

* BAM & BAI files
* FASTQ files
* SRA accessions
    <p> </p>
</body>
</html>"""


f = open('multiqc_report.html','wb')
f.write(html_template)
f.close()