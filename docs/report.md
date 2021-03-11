# Report

The report is divided into 3 sections: Summary Table, Hits Heatmap, and Hits table. An example run is available [here](https://cloudos.lifebit.ai/public/jobs/6000598e57851f0112a88fe6)

The hits are classified as such:

- **Perfect hit:**  Perfect matches to the curated reference sequences and mutations in the CARD.
- **Strict hit:** Detects previously unknown variants of known AMR genes, including secondary screen for key mutations, using detection models with CARD's curated similarity cut-offs to ensure the detected variant is likely a functional AMR gene.
- **Loose hit:** Works outside of the detection model cut-offs to provide detection of new, emergent threats and more distant homologs of AMR genes, but will also catalog homologous sequences and spurious partial hits that may not have a role in AMR.

## Summary Table

Lists the number of Perfect, Strict and Loose hits for each sample.  

## Hits Heatmap

Heatmap with RGI results for perfect and strict results. Hits are colour coded by the type of hit, ranging from blue (perfect), dark gray (strict) and white (no hit). Sample name, symbol and hit type are available as hoover information.

## Hits Table

Exhaustive table containing information on gene symbol, gene family, drug class and resistance mechanism for all the hits.
