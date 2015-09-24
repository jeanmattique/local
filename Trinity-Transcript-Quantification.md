#Trinity Transcript Quantification

There are now several methods available for estimating transcript abundance in a genome-free manner, and these include alignment-based methods (aligning reads to the transcript assembly) and alignment-free methods (typically examining k-mer abundances in the reads and in the resulting assemblies).

In Trinity, we provide direct support for running the alignment-based quantification methods [RSEM](http://deweylab.biostat.wisc.edu/rsem/) and [eXpress](http://bio.math.berkeley.edu/eXpress/), as well as the ultra-fast alignment-free method [kallisto](http://pachterlab.github.io/kallisto/).

The Trinity software does not come pre-packaged with any of these software tools, so be sure to download and install any that you wish to use. The tools should be available via your PATH setting (so, typing 'which kallisto' or 'which express' on the linux command line returns the path to where the tool is installed on your system).

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

[NOTE]
It is useful to first run 'align_and_estimate_abundance.pl' to only prep your reference database for alignment, using '--prep_reference', and then subsequently running it on each of your sets of reads in parallel to obtain sample-specific abundance estimates.

[NOTE]
If you quality-trimmed your reads using the --trimmomatic parameter in Trinity, you should consider using the corresponding quality-trimmed reads for the abundance estimation process outlined here.


## Alignment-based methods
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

