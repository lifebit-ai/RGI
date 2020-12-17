#!/usr/bin/env python3

SUMMARY_TABLE = "$summary_html"
HITS_TABLE = "$html_hit_table"
HEATMAP = "$heatmap"

html_template = """
<!DOCTYPE html>
<html>
<title>RGI-nf Report</title>
<head><meta charset="utf-8" />
</head>
<body>

<div>
<h1 id="RGI-nf-Report" style="text-align: center;"><strong>RGI-nf Report</strong></h1>
<p style="text-align: justify;"><strong>RGI-nf</strong> integrates <a href="https://github.com/arpcard/rgi">The Resistance Gene Identifier (RGI)</a> package to to predict resistome(s) from nucleotide data based on homology and SNP models. The application uses reference data from the <a href="https://card.mcmaster.ca/">Comprehensive Antibiotic Resistance Database (CARD)</a>.</p>
</div>

<div>
<ul>
<li id="Results-Summary">
<h2>Results Summary</h2>
</li>
</ul>
<p style="text-align: left;"><strong>Table 1:</strong> Hit count per sample.</p>
<p style="text-align: left;">{0}</p>
</div>

<div>
<ul>
<li id="Results-Summary">
<h2>AMR Genes</h2>
</li>
</ul>
{1}
<p style="text-align: center;"><strong>Figure 1:</strong> Heatmap of AMR genes detected per sample. Bright blue represents a perfect hit, gray represents a strict hit, white represents no hit.</p>
<p>&nbsp;</p>
</div>

<div>
<p>&nbsp;</p>
<p style="text-align: left;"><strong>&nbsp;Table 2:</strong> Drug class, resistance mechanism and gene family for the&nbsp;AMR genes detected for all samples.</p>
<p style="text-align: left;">{2}</p>
</div>

</body>
</html>
"""

def main(summary_file, hits_table, heatmap):

    # TODO - load RGI Logo?

    # Results Summary
    with open(summary_file) as summary_fh:
        summary_table_html = summary_fh.read()

    # Hits table
    with open(hits_table) as hits_table_fh:
        hits_table_html = hits_table_fh.read()

    # Heatmap
    import base64
    data_uri = base64.b64encode(open(heatmap, 'rb').read()).decode('utf-8')
    img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)

    f = open('multiqc_report.html','w')
    f.write(html_template.format(summary_table_html, img_tag, hits_table_html))
    f.close()

if __name__ == "__main__":
    main(SUMMARY_TABLE, HITS_TABLE, HEATMAP)
