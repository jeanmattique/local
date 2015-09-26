#Trinity: Frequently Asked Questions

- [There are too many transcripts! What do I do?](#ques_why_so_many_transcripts)
- [What computing resources are required?](#ques_comp_resources_required)
- [How long should this take?](#ques_how_long)
- [How do I use reads I downloaded from SRA?](#ques_sra_fq_conversion)
- [How can I run this in parallel on a computing grid?](#ques_computing_grid)
- [How do I identify the specific reads that were incorporated into the transcript assemblies?](#ques_reads_in_assembly)
- [How do I combine multiple libraries in a single Trinity run? Or how do I combine paired and single reads?](#ques_mult_seq_libraries)
- [Trinity process died due to 'std::bad_alloc'](#ques_bad_alloc)
- [Butterfly fails with java Error: Cannot create GC thread. Out of system resources.](#ques_butterfly_GC_thread_fail)
- [Can I combine strand-specific and non-strand-specific reads?](#ques_combine_SS_w_DS_reads)


<a name="ques_why_so_many_transcripts"></a>
##There are too many transcripts!  What do I do?

Welcome to RNA-Seq de novo assembly!  :)  There are several aspects of de novo assembled transcripts to be aware of:

-  *Lots* of transcripts is the rule rather than the exception.  

-  Most of the transcripts are very lowly expressed, and the deeper you sequence and the more complex your genome, the larger the number of lowly expressed transcripts you will be able to assemble.  Biological relevance of the lowly expressed transcripts could be questionable - some are bound to be very relevant.

-  There's really no good reason to filter them out.  They can be 'passengers' throughout all of your data analyses, and if any of them are important, they'll ideally surface in the relevant study.   You can put them all through Trinotate.github.io for annotation/analysis, and you can put them through DE studies just fine (those with insufficent reads will get naturally filtered out during DE studies).  If the read counts are few or lacking, they simply won't surface as significant DE entries, but if there's protein homology or other interesting features, you'll continue to capture this info.

-  If you want to guestimate 'how many expressed genes/transcripts do I really have?', then examine link:abundance_estimation.html#how_many_expr[counts of transcripts or genes vs. min expression thresholds]. This will enable you to count the number that reflect the majority of the reads, excluding counting those that have little read support.  The entries with >= ~1 fpkm or tpm tend to be heavily enriched for transcripts that correspond to what we typically think of as 'genes' in the pre-RNA-Seq era, and what typically are awarded annotation status in genome annotation. 

-  If for some reason you do want to filter out those transcripts that have little read support, you can link:abundance_estimation.html#filtering_transcripts[filter them based on the RSEM abundance estimates].

-  If you're concerned about having too many isoforms for a given gene, you can use the filtering method above to remove isoforms that are weakly expressed compared to more dominant isoforms.  Also, you can simply perform DE analysis at the 'gene' level in addition to the 'isoform' level, so as to increase your power for DE analysis.  Other methods such as [Corset](http://genomebiology.com/2014/15/7/410) can also be used to regroup relevant transcripts into 'gene' clusters.


<a name="ques_comp_resources_required"></a>
##What computing resources are required?

Ideally, you will have access to a large-memory server, roughly having ~1G of RAM per 1M reads to be assembled.  The memory usage mostly depends on the complexity of the RNA-Seq data set, specifically on the number of unique k-mers.  If you do not have access to a high-memory server, link:index.html#RunElsewhere[other freely available options are available].

<a name="ques_how_long"></a>
##How long should this take?

Trinity run-time depends on a number of factors including the number of reads to be assembled and the complexity of the transcript graphs.  The assembly from start to finish can take anywhere from ~1/2 hour to 1 hour per million reads (your mileage may vary).  If you have a large data set, be sure to use the --normalize_reads parameter to greatly improve run times.



<a name="ques_sra_fq_conversion"></a>
##How do I use reads I downloaded from SRA?

RNA-Seq data downloaded from SRA tends to exist in a .sra file that needs to be converted to fastq file format.  This can be done using the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software) like so:

      SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra

<a name="ques_computing_grid"></a>
##How can I run this in parallel on a computing grid?

The initial Inchworm and Chrysalis steps currently need to be run on a single server, however parts of Chrysalis and all of Butterfly can be run as naive parallel processes on a computing grid. Both SGE and LSF are currently supported. Configuring Trinity for leveraging these grid-computing platforms is described [here](Installing-Trinity#grid_conf.


<a name="ques_reads_in_assembly"></a>
##How do I identify the specific reads that were incorporated into the transcript assemblies?

Currently, the mappings of reads to transcripts are not reported.  To obtain this information, we recommend realigning the reads to the assembled transcripts using http://bowtie-bio.sourceforge.net/index.shtml[Bowtie]. Capturing both properly paired read alignments and single unpaired alignments can be captured and quantified as described link:analysis/abundance_estimation.html#detailed_assessment[here].


<a name="ques_mult_seq_libraries"></a>
##How do I combine multiple libraries in a single Trinity run? Or, how do I combine paired and single reads?

If you have RNA-Seq data from multiple libraries and you want to run them all through Trinity in a single pass, simply combine all your left.fq files into one left.fq file, and combine all right.fq files into one right.fq file. Then run Trinity using these separately concatenated left and right input files.  

If you have additional singletons, add them to the .fq file that they correspond to based on the sequencing method used (if they're equivalent to the left.fq entries, add them there, etc).

There is no good way to combine strand-specific data with non-strand-specific data, unless you decide to treat the entire data set as non-strand-specific.


<a name="ques_bad_alloc"></a>
##Trinity process died due to 'std::bad_alloc'

This is an indicator that the process ran out of available RAM. If you have more RAM resources to make available to Trinity, then simply rerun your original Trinity command with the altered resources allocated and it should resume approximately where it left off.  

If you are resource limited, please consider running Trinity link:index.html#RunElsewhere[here].  If you want to continue to try to run Trinity given your available resources, you can reduce the total RAM requirements by running Trinity with parameter '--min_kmer_cov 2'. Although the assembly should still be of high quality and require less RAM, lowly expressed transcripts may be more highly fragmented in the assembly.


<a name="ques_butterfly_GC_thread_fail"></a>
##Butterfly fails with java Error: Cannot create GC thread. Out of system resources.

There are a couple reasons why this error message might creep up.

1.  *all memory has been consumed on the machine*.  Each butterfly process wants to reserve 10G of maximum heap space.  If there's less than 10G of free memory on the machine per butterfly (--CPU setting), then java may not be able to initialize (depends on your OS configuration).  Try reducing the --CPU setting and rerunning your Trinity command. It should resume where it left off.

2.  *NUMA architecture*:  one of our users found that the java invocation required: -XX:ParallelGCThreads=<Numerical Thread Count>, otherwise it would try to use too many threads.

<a name="ques_combine_SS_w_DS_reads"></a>
##Can I combine strand-specific and non-strand-specific reads? 

You can do so, but you wouldn't be able to benefit as from the strand-specificity, since you'll need to run Trinity in non-strand-specific mode.
