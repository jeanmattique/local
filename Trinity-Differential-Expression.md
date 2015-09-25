# Differential Expression Analysis Using a Trinity Assembly

Our current system for identifying differentially expressed transcripts relies on using the EdgeR Bioconductor package. We have a protocol and scripts described below for identifying differentially expressed transcripts and clustering transcripts according to expression profiles. This process is somewhat interactive, and described are automated approaches as well as manual approaches to refining gene clusters and examining their corresponding expression patterns.

We recommend generating a single Trinity assembly based on combining all reads across all samples as inputs.  Then, reads are separately aligned back to the single Trinity assembly for downstream analyses of differential expression, according to our link:abundance_estimation.html[abundance estimation protocol].   If you decide to assemble each sample separately, then you'll likely have difficulty comparing the results across the different samples due to differences in assembled transcript lengths and contiguity.

Before attempting to analyze differential expression, you should have already [estimated transcript abundance and generated an RNA-Seq counts matrix containing RNA-Seq fragment counts](Trinity-Transcript-Quantification) for each of your transcripts (or genes) across each biological replicate for each sample (experiment, condition, tissue, etc.).


