#Trinity Transcript Quantification

There are now several methods available for estimating transcript abundance in a genome-free manner, and these include alignment-based methods (aligning reads to the transcript assembly) and alignment-free methods (typically examining k-mer abundances in the reads and in the resulting assemblies).

In Trinity, we provide direct support for running the alignment-based quantification methods [RSEM](http://deweylab.biostat.wisc.edu/rsem/) and [eXpress](http://bio.math.berkeley.edu/eXpress/), as well as the ultra-fast alignment-free method [kallisto](http://pachterlab.github.io/kallisto/).

The Trinity software does not come pre-packaged with any of these software tools, so be sure to download and install any that you wish to use. The tools should be available via your PATH setting (so, typing 'which kallisto' or 'which express' on the linux command line returns the path to where the tool is installed on your system).

>If you have multiple RNA-Seq data sets that you want to compare (eg. different tissues sampled from a single organims), be sure to **generate a single Trinity assembly** and to then run the abundance estimation separately for each of your samples.

## Estimating Transcript Abundance

The Trinity toolkit comes with a script to facilitate running your choice of the above tools to quantitate transcript abundance:

     % $TRINITY_HOME/util/align_and_estimate_abundance.pl 

     #########################################################################
    #
    #  --transcripts <string>           transcript fasta file
    #  --seqType <string>               fq|fa
    # 
    #  If Paired-end: 
    #
    #  --left <string>
    #  --right <string>
    #  
    #    or Single-end:
    #
    #  --single <string>
    #
    #  --est_method <string>           abundance estimation method.
    #                                        alignment_based:  RSEM|eXpress       
    #                                        alignment_free: kallisto
    #  
    # --output_dir <string>            write all files to output directory
    #  
    #
    #  if alignment_based est_method:
    #       --aln_method <string>            bowtie|bowtie2|(path to bam file) alignment method.  (note: RSEM requires bowtie)
    #                                       (if you already have a bam file, you can use it here instead of rerunning bowtie)
    #
    # Optional:
    # 
    # --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
    #
    # --thread_count                   number of threads to use (default = 4)
    #
    # --max_ins_size <int>             maximum insert size (bowtie -X parameter, default: 800)
    #
    # --debug                          retain intermediate files
    #
    #
    #  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
    #     or  
    #  --trinity_mode                   Setting --trinity_mode will automatically generate the gene_trans_map and use it.
    #
    #
    #  --prep_reference                 prep reference set for eXpress (builds bowtie index, etc)
    #
    #  --output_prefix <string>         prefix for output files.  Defaults to --est_method setting.
    #
    #
    #  if alignment_based method:
    #        --coordsort_bam                  provide coord-sorted bam in addition to the default (unsorted) bam.
    #
    #  --show_full_usage_info           provide more detailed usage info for customizing the alignment or abundance estimation parameters.
    #
    #############################
    #  RSEM opts:
    #  --fragment_length <int>         optionally specify fragment length (not seq length, but frag size ie. 300) for SE reads.
    #
    #  --include_rsem_bam              provide the RSEM enhanced bam file including posterior probabilities of read assignments.
    #
    #########################################################################
    #
    #  Example usage:
    #
    #   ## Just prepare the reference for alignment and abundance estimation
    #
    #    ./util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
    #
    #   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
    #
    #    ./util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode 
    #
    ##  ## prep the reference and run the alignment/estimation
    #
    #    ./util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
    #
    #########################################################################


If you have strand-specific data, be sure to include the '--SS_lib_type' parameter.

>It is useful to first run 'align_and_estimate_abundance.pl' to only prep your reference database for alignment, using '--prep_reference', and then subsequently running it on each of your sets of reads in parallel to obtain sample-specific abundance estimates.

>If you quality-trimmed your reads using the --trimmomatic parameter in Trinity, you should consider using the corresponding quality-trimmed reads for the abundance estimation process outlined here.


Any of the abundance estimation methods will provide transcript-level estimates of the count of RNA-Seq fragments that were derived from each transcript, in addition to a normalized measure of transcript expression that takes into account the transcript length, the number of reads mapped to the transcript, and the the total number of reads that mapped to any transcript. Normalized expression metrics may be reported as 'fragments per kilobase transcript length per million fragments mapped' (FPKM) or 'transcripts per million transcripts' (TPM).  The TPM metric is generally preferred to FPKM, given the property that all values will always sum up to 1 million (FPKM values will tend to not sum up to the same value across samples). For info on TPM vs. FPKM, see [Wagner et al. 2012. Theory Biosci](http://lynchlab.uchicago.edu/publications/Wagner,%20Kin,%20and%20Lynch%20(2012).pdf), and [Li and Dewey, BMC Bioinf. 2011](http://www.biomedcentral.com/1471-2105/12/323).

The estimated fragment counts are generally needed by many differential expression analysis tools that use count-based statistical models, and the normalized expression values (FPKM or TPM) are used almost everywhere else, such as plotting in heatmaps.


## Alignment-based abundance estimation methods
The alignment step generates the file 'bowtie.bam', which is then fed directly into either RSEM or eXpress.  Note, if parameter '--coordsort_bam ' is set, the process also generates a 'bowtie.csorted.bam' file, which is a coordinate-sorted bam file that can be used for visualization using IGV.

### RSEM output
The RSEM computation generates two primary output files containing the abundance estimation information:

  RSEM.isoforms.results  : EM read counts per Trinity transcript
  RSEM.genes.results     : EM read counts on a per-Trinity-gene, 'gene' used loosely here.


The output for the isoforms file looks like so:

|transcript_id|   gene_id| length|  effective_length|        expected_count|  TPM|     FPKM|    IsoPct|
|-------------:|---------:|-----:|-----------------:|---------------------:|-----:|--------:|--------:|
|TRINITY_DN100_c0_g1_i1|  TRINITY_DN100_c0_g1|     443|     181.37|  4.84|    1670.06| 12311.85|        100.00|
|TRINITY_DN101_c0_g1_i1|  TRINITY_DN101_c0_g1|     251|     19.37|   1.00|    3231.22| 23820.87|        100.00|
|TRINITY_DN103_c0_g1_i1|  TRINITY_DN103_c0_g1|     1219|    957.37|  0.00|    0.00|    0.00|    100.00|
|TRINITY_DN103_c0_g2_i1|  TRINITY_DN103_c0_g2|     414|     152.41|  0.00|    0.00|    0.00|    100.00|
|TRINITY_DN104_c0_g1_i1|  TRINITY_DN104_c0_g1|     408|     146.44|  0.00|    0.00|    0.00|    0.00|
|TRINITY_DN106_c0_g1_i1|  TRINITY_DN106_c0_g1|     1111|    849.37|  1.00|    73.70|   543.30|  100.00|
|TRINITY_DN106_c1_g1_i1|  TRINITY_DN106_c1_g1|     339|     81.68|   0.00|    0.00|    0.00|    0.00|


### eXpress output

The eXpress runner also generates two files:

     results.xprs  : transcript abundance estimates (generated by eXpress)
     results.xprs.genes : gene abundance estimates (provided here by summing up transcript values per gene)

Each includes the estimated counts, FPKM, and TPM measures, in addition to a large number of other metrics - see the [eXpress documentation](http://bio.math.berkeley.edu/eXpress/manual.html) for details.

## Alignment-free abundance estimation methods

### kallisto

The kallisto runner similarly generates two files:

     abundance.tsv  : transcript abundance estimates (generated by kallisto)
     abundance.tsv.genes : gene abundance estimates (provided here by summing up transcript values per gene)

The format of the output is short and sweet, providing the key essential metrics:

|target_id|       length|  eff_length|      est_counts|      tpm|
|:---------|-----------:|-----------:|---------------:|--------:|
|TRINITY_DN10_c0_g1_i1|   334|     100.489| 13|      4186.62|
|TRINITY_DN11_c0_g1_i1|   319|     87.9968| 0|       0|
|TRINITY_DN12_c0_g1_i1|   244|     38.2208| 2|       1693.43|
|TRINITY_DN17_c0_g1_i1|   229|     30.2382| 5|       5351.21|
|TRINITY_DN18_c0_g1_i1|   633|     384.493| 19|      1599.2|
|TRINITY_DN18_c1_g1_i1|   289|     65.795|  1|       491.864|
|TRINITY_DN19_c0_g1_i1|   283|     61.0618| 10|      5299.91|


## Build Transcript and Gene Expression Matrices

Using the transcript and gene-level abundance estimates for each of your samples, construct a matrix of counts and a matrix of normalized expression values using the following script:

    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl 

    ############################################################
    #
    # Usage:  $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method <method>  sample1.results sample2.results ...
    # Required:
    #
    #  --est_method <string>           RSEM|eXpress|kallisto  (needs to know what format to expect)
    #
    # Options:
    #
    #  --cross_sample_norm <string>         TMM|UpperQuartile|none   (default: TMM)
    #
    #  --name_sample_by_basedir             name sample column by dirname instead of filename
    #      --basedir_index <int>            default(-2)
    #
    #  --out_prefix <string>                default: 'matrix'
    #
    ############################################################



For example, suppose you have two samples (sampleA and sampleB), and you ran kallisto to estimate transcript abundances, you might generate matrices like so:

    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto \
        --out_prefix trans_counts \
        --name_sample_by_basedir \
         sampleA/abundance.tsv \
         sampleB/abundance.tsv 

which would generate the following three files:

      trans_counts.counts.matrix  : the estimated RNA-Seq fragment counts (raw counts)
      trans_counts.TPM.not_cross_norm  : a matrix of TPM expression values (not cross-sample normalized)
      trans_counts.TMM.EXPR.matrix : a matrix of TMM-normalized expression values

The 'counts.matrix' file is used for downstream analyses of differential expression.  The TMM.EXPR.matrix file is used as the gene expression matrix in most other analyses.  For information on the importance of TMM (or cross-sample normalization in general), see [Robinson & Oshlack, Genome Biology 2010](http://www.genomebiology.com/2010/11/3/R25) and [Dillies et al., Brief Bioinf, 2012](http://bib.oxfordjournals.org/content/14/6/671.long).
