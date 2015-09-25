# Differential Expression Analysis Using a Trinity Assembly

Our current system for identifying differentially expressed transcripts relies on using the EdgeR Bioconductor package. We have a protocol and scripts described below for identifying differentially expressed transcripts and clustering transcripts according to expression profiles. This process is somewhat interactive, and described are automated approaches as well as manual approaches to refining gene clusters and examining their corresponding expression patterns.

We recommend generating a single Trinity assembly based on combining all reads across all samples as inputs.  Then, reads are separately aligned back to the single Trinity assembly for downstream analyses of differential expression, according to our link:abundance_estimation.html[abundance estimation protocol].   If you decide to assemble each sample separately, then you'll likely have difficulty comparing the results across the different samples due to differences in assembled transcript lengths and contiguity.

Before attempting to analyze differential expression, you should have already [estimated transcript abundance and generated an RNA-Seq counts matrix containing RNA-Seq fragment counts](Trinity-Transcript-Quantification) for each of your transcripts (or genes) across each biological replicate for each sample (experiment, condition, tissue, etc.).

>If you have biological replicates, be sure to align each replicate set of reads and estimate abundance values for the sample independently, and targeting the single same targeted Trinity assembly.

## Running Differential Expression Analysis

Trinity provides support for several differential expression analysis tools, currently including the following R packages:

* edgeR : <http://bioconductor.org/packages/release/bioc/html/edgeR.html>
* DESeq2: <http://bioconductor.org/packages/release/bioc/html/DESeq2.html>
* limma/voom: <http://bioconductor.org/packages/release/bioc/html/limma.html>
* ROTS: <http://www.btk.fi/research/research-groups/elo/software/rots/>

Be sure to have [R](https://www.r-project.org/) installed in addition to the above software package that you want to use for DE detection. 

In addition, you'll need the following R packages installed: ctc, Biobase, gplots, and ape.  These can be installed like so (along with the subset of Bioconductor packages for DE software above)

     % R
     > source("http://bioconductor.org/biocLite.R")
     > biocLite('edgeR')
     > biocLite('limma')
     > biocLite('DESeq2')
     > biocLite('ctc')
     > biocLite('Biobase')
     > install.packages('gplots')
     > install.packages('ape')


Differentially expressed transcripts or genes are identified by running the script below, which will perform pairwise comparisons among each of your sample types. If you have biological replicates for each sample, you should indicate this as well (described further below).  To analyze transcripts, use the 'transcripts.counts.matrix' file. To perform an analysis at the 'gene' level, use the 'genes.counts.matrix'. Again, Trinity Components are used as a proxy for 'gene' level studies.

     $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl 

     #################################################################################################
     #
     #  Required:
     #
     #  --matrix|m <string>               matrix of raw read counts (not normalized!)
     #
     #  --method <string>               edgeR|DESeq2|voom|ROTS
     #                                     note: you should have biological replicates.
     #                                           edgeR will support having no bio replicates with
     #                                           a fixed dispersion setting. 
     #
     #  Optional:
     #
     #  --samples_file|s <string>         tab-delimited text file indicating biological replicate relationships.
     #                                   ex.
     #                                        cond_A    cond_A_rep1
     #                                        cond_A    cond_A_rep2
     #                                        cond_B    cond_B_rep1
     #                                        cond_B    cond_B_rep2
     #
     #
     #  General options:
     #
     #  --min_rowSum_counts <int>       default: 2  (only those rows of matrix meeting requirement will be tested)
     #
     #  --output|o                      name of directory to place outputs (default: $method.$pid.dir)
     #
     #  --reference_sample <string>     name of a sample to which all other samples should be compared.
     #                                   (default is doing all pairwise-comparisons among samples)
     #
     #  --contrasts <string>            file (tab-delimited) containing the pairs of sample comparisons to perform.
     #                                  ex. 
     #                                       cond_A    cond_B
     #                                       cond_Y    cond_Z
     #
     #
     ###############################################################################################
     #
     #  ## EdgeR-related parameters
     #  ## (no biological replicates)
     #
     #  --dispersion <float>            edgeR dispersion value (Read edgeR manual to guide your value choice)
     #
     #  ## ROTS parameters
     #  --ROTS_B <int>                   : number of bootstraps and permutation resampling (default: 500)
     #  --ROTS_K <int>                   : largest top genes size (default: 5000)
     # 
     #
     ###############################################################################################
     #
     #   Documentation and manuals for various DE methods.  Please read for more advanced and more
     #   fine-tuned DE analysis than provided by this helper script.
     #
     #  edgeR:       http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
     #  DESeq2:      http://bioconductor.org/packages/release/bioc/html/DESeq2.html    
     #  voom/limma:  http://bioconductor.org/packages/release/bioc/html/limma.html
     #  ROTS:        http://www.btk.fi/research/research-groups/elo/software/rots/
     #
     ###############################################################################################



### Identifying DE Features: No Biological Replicates (Proceed with Caution)

It's very important to have biological replicates to power DE detection and reduce false positive predictions. If you do not have biological replicates, edgeR will allow you to perform DE analysis if you manually set the --dispersion parameter. Values for the dispersion parameter must be chosen carefully, and you might begin by exploring values between 0.1 and 0.4. Please visit the [edgeR manual](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) for further guidance on this matter.

      $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
          --matrix counts.matrix --method edgeR --dispersion 0.1


### Identifying DE features: With biological replicates (PREFERRED)

Be sure to have a tab-delimited 'samples_described.txt' file that describes the relationship between samples and replicates.  For example:

    conditionA   condA-rep1
    conditionA   condA-rep2
  
    conditionB   condB-rep1
    conditionB   condB-rep2
  
    conditionC   condC-rep1
    conditionC   condC-rep2


where condA-rep1, condA-rep2, condB-rep1, etc..., are all column names in the 'counts.matrix' generated earlier (see top of page). Your sample names that group the replicates are user-defined here.


Any of the available methods support analyses containing biological replicates.  Here, for example, we again choose voom within the limma package.

      $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
       --matrix counts.matrix --method voom --samples_file samples_described.txt 
    

>A full example of the edgeR pipeline involving combining reads from multiple samples, assembling them using Trinity, separately aligning reads back to the trintiy assemblies, abundance estimation using RSEM, and differential expression analysis using edgeR is provided at: $TRINITY_HOME/sample_data/test_full_edgeR_pipeline
