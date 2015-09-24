# RNA-Seq De novo Assembly Using Trinity

![TrinityCompositeLogo](https://raw.githubusercontent.com/wiki/trinityrnaseq/trinityrnaseq/images/TrinityCompositeLogo.png)

## Quick Guide for the Impatient

Trinity assembles transcript sequences from Illumina RNA-Seq data.

Download Trinity [here](https://github.com/trinityrnaseq/trinityrnaseq/releases).

Build Trinity by typing 'make' in the base installation directory.

Assemble RNA-Seq data like so:

   Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G 

Find assembled transcripts as:  'trinity_out_dir/Trinity.fasta'


== Intro to Trinity ==

Trinity, developed at the http://www.broadinstitute.org[Broad Institute] and the http://www.cs.huji.ac.il[Hebrew University of Jerusalem], represents a novel method for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data. Trinity combines three independent software modules: Inchworm, Chrysalis, and Butterfly, applied sequentially to process large volumes of RNA-seq reads. Trinity partitions the sequence data into many individual de Bruijn graphs, each representing the transcriptional complexity at at a given gene or locus, and then processes each graph independently to extract full-length splicing isoforms and to tease apart transcripts derived from paralogous genes.  Briefly, the process works like so:

- *Inchworm* assembles the RNA-seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.

- *Chrysalis* clusters the Inchworm contigs into clusters and constructs complete de Bruijn graphs for each cluster.  Each cluster represents the full transcriptonal complexity for a given gene (or sets of genes that share sequences in common).  Chrysalis then partitions the full read set among these disjoint graphs.

- *Butterfly* then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.

Trinity was published in http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571712/[Nature Biotechnology].  Our protocol for transcriptome assembly and downstream analysis is now published in http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/[Nature Protocols]

The Trinity software package can be downloaded https://github.com/trinityrnaseq/trinityrnaseq/releases[here on GitHub]. Legacy versions (pre-2015) are still available at http://sourceforge.net/projects/trinityrnaseq/files/PREV_CONTENTS/previous_releases/[Sourceforge].

http://trinityrnaseq.github.io/performance/[Runtime and transcript reconstruction performance stats] are available for current and previous releases.

http://www.broadinstitute.org/partnerships/education/broade/trinity-screencast[Screencast videos] are available to introduce you to Trinity and its various components. Also, hands-on tutorials for Trinity and Tuxedo are available as part of our link:workshop/rnaseq_workshop.html[RNA-Seq Workshop].




== Table of Contents ==

* <<installation, Installing Trinity>>
* <<running_trinity, Running Trinity>>
** <<typical_usage, Typical Trinity Command Line>>
** <<typical_options, Options to consider when running Trinity>>
*** <<trimmomatic, Quality trimming using Trimmomatic>>
*** <<insilinorm, In silico Normalization of Reads (for assembly of (hundreds of millions to billions of reads)>>
*** <<jaccard_clip, Minimizing falsely fused transcripts from gene-dense genomes>>
*** <<genome_guided, Genome-guided Trinity>>
*** <<genome_annotation, Comprehensive transcriptome-based genome annotation using Trinity and PASA>>
** <<trinity_output, Output of Trinity>>
** <<compute_requirements, Hardware and Configuration Requirements>>
** <<monitoring_trinity, Monitoring the Progress of Trinity>>
** <<sample_data, Running Trinity on Sample Data>>
* <<Downstream_analyses, Post-assembly Transcriptome Analysis>>
** link:analysis/abundance_estimation.html[Abundance Estimation using RSEM or eXpress, and Visualization using IGV]
** link:analysis/diff_expression_analysis.html[Differential Expression Analysis using Bioconductor]
** http://transdecoder.github.io[Protein-coding Region Identification Using TransDecoder]
** http://trinotate.github.io[Functional Annotation Using Trinotate]
** link:analysis/run_GOseq.html[Gene Ontology functional category enrichment using GOseq and Trinotate]
** link:analysis/full_length_transcript_analysis.html[Full-length Transcript Analysis]
* link:advanced_trinity_guide.html[Advanced Guide to Trinity]
* link:trinity_faq.html[Frequently Asked Questions]
* <<trinity_tidbits, Trinity Tidbits>>
* <<contact_us, Contact Us>>
* <<referencing_trinity, Referencing Trinity>>


[[installation]]
== Installing Trinity ==

=== Local Installation of Trinity on a High-memory Linux Server ===

After https://github.com/trinityrnaseq/trinityrnaseq/releases[downloading] the sofware to a Linux server, simply type 
   
   make 

in the base installation directory.  This should build Inchworm and Chrysalis, both written in C++.  Butterfly should not require any special compilation, as its written in Java and already provided as portable precompiled software, but *Java-1.7* is required.

Afterwards, you may want to build the additional plugin components that provide support for downstream analyses (such as abundance estimation using RSEM), in which case you would then type:

   make plugins

[NOTE]
==================
If you encounter any errors in building the RSEM software, simply

   cd trinity-plugins/tmp.rsem

   make

and assuming that succeeds, you can then cd back to the main Trinity installation directory and retype 'make plugins' to continue the remaining build.
==================

Additional tools required for running Trinity include:

- http://bowtie-bio.sourceforge.net/index.shtml[bowtie-1]

Trinity has been tested and is supported on Linux.


[[Computing_Grid]]
== Adapting Trinity to a computing grid for parallel processing of naively parallel steps ==

[NOTE]
Trinity supports LSF, SGE, SLURM, and PBS.

Trinity has many parallel-components, all of which can benefit from having multiple CPUs on a single server, but there are also cases such as in Chrysalis and Butterfly where tens of thousands to hundreds of thousands of commands can be executed in parallel, each having independent inputs and outputs.  These naively-parallel commands can be most efficiently computed in the context of a compute farm, submitting each of the commands (or batches of them) to individual nodes on the computing grid.  There are several different computing grid job management systems that are in common use, such as SGE or LSF.

Trinity currently supports both SGE and LSF.  To leverage either, simply run 'Trinity --grid_conf your_conf_file.txt', where your_conf_file.txt is a very simple configuration file that indicates parameters for the grid job submission. For example, at the Broad and using LSF, a configuration file might contain the following:

 #-------------------------------------------------------------------------------------------
 # grid type: 
 grid=LSF
 
 # template for a grid submission
 cmd=bsub -q regevlab -R "rusage[mem=10]"
 # note -e error.file -o out.file are set internally, so dont set them in the above cmd. 
 
 # uses the LSF feature to pre-exec and check that the file system is mounted before executing.
 # this helps when you have some misbehaving grid nodes that lost certain file mounts.
 mount_test=T
 
 ##########################################################################################
 # settings below configure the Trinity job submission system, not tied to the grid itself.
 ##########################################################################################
 
 # number of grid submissions to be maintained at steady state by the Trinity submission system 
 max_nodes=500
 
 # number of commands that are batched into a single grid submission job.
 cmds_per_node=100

 #--------------------------------------------------------------------------------------------


where the above indicates that LSF is the grid type (either LSF or SGE are supported), the queue to submit to is our 'regevlab' named queue, and memory is set to 10 gigabytes. Up to 500 jobs will be submitted at any given time (throttled by the Trinity-included job management system), and the jobs are batched at 10 commands per submission (so, for example, 10 butterfly jobs will be submitted as a single grid job, each being executed serially).

For SGE, at the Broad Institute, we might specify a configuration:

 #--------------------------------------------------------------------------------------------
 # grid type: 
 grid=SGE
 # template for a grid submission
 cmd=qsub -V -cwd
 # number of grid submissions to be maintained at steady state by the Trinity submission system 
 max_nodes=500
 # number of commands that are batched into a single grid submission job.
 cmds_per_node=1
 #--------------------------------------------------------------------------------------------

where, SGE is indicated as the grid type.  We don't need to specify a queue name, apparently, as it gets submitted to the default queue, and the default memory allocation is sufficient. The project_code can also be left blank unless your SGE configuration requires it.  The maximum number of nodes to throttle the jobs at (500) and the number of commands executed in a single grid job (10) is the same as what we show above for our LSF configuration.

Likewise, for SLURM, we have:

 #---------------------------------------------------------------------------------------------
 # grid type: 
 grid=SLURM
 # template for a grid submission
 cmd=sbatch -p queue_name --mem=10000 --time=02:00:00 
 # number of grid submissions to be maintained at steady state by the Trinity submission system 
 max_nodes=4000
 # number of commands that are batched into a single grid submission job.
 cmds_per_node=20
 #----------------------------------------------------------------------------------------------


Example configuration files are provided under $TRINITY_HOME/htc_conf


[[RunElsewhere]]
=== Using a Freely Available Trinity Installation on High Performance Computing Systems ===

- Use the Trinity NCGAS Galaxy portal at https://galaxy.ncgas-trinity.indiana.edu/[https://galaxy.ncgas-trinity.indiana.edu/].

- Trinity is available on XSEDE's Blacklight server at the http://www.psc.edu/[Pittsburgh Supercomputer Center].  Information on how researchers in the USA can get a FREE account and to run Trinity on Blacklight (which has up to 16TB of RAM!) is provided http://trinity-use-on-blacklight-psc.wikispaces.com/Trinity+Usage+on+Blacklight[here]. Thanks to Phil Blood and Brian Cougar for maintaining this installation and making services available.

- http://diagcomputing.org/[The Data Intensive Acadmeic Grid (DIAG)] provides *FREE ACCESS TO ALL RESEARCHERS* high memory servers and data storage for academic research. Trinity is supported as one of the pre-installed applications. The guide for running Trinity on DIAG is http://wiki.diagcomputing.org/index.php/Trinity[here]. Thanks to Anup Mahurkar and Joshua Orvis for support.


[[running_trinity]]
== Running Trinity ==

Trinity is run via the script: 'Trinity' found in the base installation directory.

Usage info is as follows:



 ###############################################################################
 #
 #     ______  ____   ____  ____   ____  ______  __ __
 #    |      ||    \ |    ||    \ |    ||      ||  |  |
 #    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
 #    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
 #      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
 #      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
 #      |__|  |__|\_||____||__|__||____|  |__|  |____/
 #
 ###############################################################################
 #
 # Required:
 #
 #  --seqType <string>      :type of reads: ( fa, or fq )
 #
 #  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
 #                            provied in Gb of RAM, ie.  '--max_memory 10G'
 #
 #  If paired reads:
 #      --left  <string>    :left reads, one or more file names (separated by commas, not spaces)
 #      --right <string>    :right reads, one or more file names (separated by commas, not spaces)
 #
 #  Or, if unpaired reads:
 #      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
 #
 ####################################
 ##  Misc:  #########################
 #
 #  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
 #                                   if paired: RF or FR,
 #                                   if single: F or R.   (dUTP method = RF)
 #                                   See web documentation.
 #
 #  --CPU <int>                     :number of CPUs to use, default: 2
 #  --min_contig_length <int>       :minimum assembled contig length to report
 #                                   (def=200)
 #
 #  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
 #
 #  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
 #                                   (see genome-guided param section under --show_full_usage_info)
 #
 #  --jaccard_clip                  :option, set if you have paired reads and
 #                                   you expect high gene density with UTR
 #                                   overlap (use FASTQ input file format
 #                                   for reads).
 #                                   (note: jaccard_clip is an expensive
 #                                   operation, so avoid using it unless
 #                                   necessary due to finding excessive fusion
 #                                   transcripts w/o it.)
 #
 #  --trimmomatic                   :run Trimmomatic to quality trim reads
 #                                        see '--quality_trimming_params' under full usage info for tailored settings.
 #                                  
 #
 #  --normalize_reads               :run in silico normalization of reads. Defaults to max. read coverage of 50.
 #                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
 #     
 #
 #  --output <string>               :name of directory for output (will be
 #                                   created if it doesn't already exist)
 #                                   default( your current working directory: "/Users/bhaas/GITHUB/trinityrnaseq/trinity_out_dir" 
 #                                    note: must include 'trinity' in the name as a safety precaution! )
 #  
 #  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
 #
 #  --cite                          :show the Trinity literature citation
 #
 #  --version                       :reports Trinity version (BLEEDING_EDGE) and exits.
 #
 #  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
 #
 #
 ###############################################################################
 #
 #  *Note, a typical Trinity command might be:
 #
 #        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
 #
 #
 #    and for Genome-guided Trinity:
 #
 #        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
 #                --genome_guided_max_intron 10000 --CPU 6
 #
 #     see: /Users/bhaas/GITHUB/trinityrnaseq/sample_data/test_Trinity_Assembly/
 #          for sample data and 'runMe.sh' for example Trinity execution
 #
 #     For more details, visit: http://trinityrnaseq.github.io
 #
 ###############################################################################


[NOTE]
Trinity performs best with strand-specific data, in which case sense and antisense transcripts can be resolved.  For protocols on strand-specific RNA-Seq, see: http://www.ncbi.nlm.nih.gov/pubmed/21943893[Borodina T, Adjaye J, Sultan M. A strand-specific library preparation protocol for RNA sequencing. Methods Enzymol. 2011;500:79-98. PubMed PMID: 21943893].


If you have strand-specific data, specify the library type.  There are four library types:

- Paired reads:
    * *RF*: first read (/1) of fragment pair is sequenced as anti-sense (reverse(*R*)), and second read (/2) is in the sense strand (forward(*F*)); typical of the dUTP/UDG sequencing method.
    * *FR*: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the antisense strand (reverse)

- Unpaired (single) reads:
    * *F*: the single read is in the sense (forward) orientation
    * *R*: the single read is in the antisense (reverse) orientation

By setting the *--SS_lib_type* parameter to one of the above, you are indicating that the reads are strand-specific.  By default, reads are treated as not strand-specific.

image:http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/bin/nihms-537313-f0006.jpg[strand-specific library type]

Other important considerations:

- Whether you use Fastq or Fasta formatted input files, be sure to keep the reads oriented as they are reported by Illumina, if the data are strand-specific. This is because, Trinity will properly orient the sequences according to the specified library type.  If the data are not strand-specific, now worries because the reads will be parsed in both orientations.

- If you have both paired and unpaired data, and the data are NOT strand-specific, you can combine the unpaired data with the left reads of the paired fragments.  Be sure that the unpaired reads have a /1 as a suffix to the accession value similarly to the left fragment reads.  The right fragment reads should all have /2 as the accession suffix.  Then, run Trinity using the --left and --right parameters as if all the data were paired.

- If you have multiple paired-end library fragment sizes, set the '--group_pairs_distance' according to the larger insert library.  Pairings that exceed that distance will be treated as if they were unpaired by the Butterfly process.  

- by setting the '--CPU option', you are indicating the maximum number of threads to be used by processes within Trinity. Note that Inchworm alone will be internally capped at 6 threads, since performance will not improve for this step beyond that setting)


[[typical_usage]]
== Typical Trinity Command Line == 

A typical Trinity command for assembling non-strand-specific RNA-seq data would be like so, running the entire process on a single high-memory server (aim for ~1G RAM per ~1M ~76 base Illumina paired reads, but often *much* less memory is required):

Run Trinity like so:

   Trinity --seqType fq --max_memory 50G --left reads_1.fq.gz  --right reads_2.fq.gz --CPU 6

If you have multiple sets of fastq files, such as corresponding to multiple tissue types or conditions, etc., you can indicate them to Trinity like so:

   Trinity --seqType fq --max_memory 50G  --left condA_1.fq.gz,condB_1.fq.gz,condC_1.fq.gz --right condA_2.fq.gz,condB_2.fq.gz,condC_2.fq.gz --CPU 6  

Also note that fastq files can be gzip-compressed as shown above, in which case they should require a '.gz' extension.

Example data and sample pipeline are provided and described <<sample_data, here>>.

[[typical_options]]
== Options to Consider when Running Trinity ==

Trinity includes additional options to automate various aspects of RNA-Seq read processing that should be considered prior to executing the de novo assembly. This includes quality trimming of reads (using http://www.usadellab.org/cms/?page=trimmomatic[Trimmomatic]), or in silico normalization of the total reads to reduce the number of reads that are subject to de novo assembly, improving on assembly run-time.  Also, if transcripts are derived from a compact genome where overlapping UTRs are common, options are provided to mitigate the assembly of falsely end-to-end fused transcripts by analyzing the consistency of the read pairings across the length of the transcripts. These options are each detailed below.

[[trimmomatic]]
=== Quality trimming using Trimmomatic ===
To perform quality trimming of inputted fastq files, use 'Trinity --trimmomatic'.  The default settings for quality trimming are described under the full usage info for Trinity (use 'Trinity --show_full_usage_info' for complete usage info):

 ################################################################################
 #### Quality Trimming Options ####  
 # 
 #  --quality_trimming_params <string>   defaults to: "LEADING:5 TRAILING:5 MINLEN:36"
 #
 ################################################################################

The various options that are available for the Trimmomatic software are described on the http://www.usadellab.org/cms/?page=trimmomatic[Trimmomatic software website].  The Trimmomatic software is bundled as a trinity plugin for convenience.


[[insilinorm]]
== Assembling Large RNA-Seq Data Sets (hundreds of millions to billions of reads) ==

If you have especially large RNA-Seq data sets involving many hundreds of millions of reads to billions of reads, consider performing an in silico normalization of the full data set using 'Trinity --normalize_reads'.  The default normalization process should work well for most data sets. If you prefer to manually set normalization-related parameters, you can find the options under the full Trinity usage info:

 ################################################################################
 ####  In silico Read Normalization Options ###
 #
 #  --normalize_max_read_cov <int>       defaults to 50
 #  --normalize_by_read_set              run normalization separate for each pair of fastq files,
 #                                       then one final normalization that combines the individual normalized reads.
 #                                       Consider using this if RAM limitations are a consideration.
 #
 ################################################################################


If you are interested in running the normalization utility outside of Trinity, you can run it directly as described link:trinity_insilico_normalization.html[here].  

[[jaccard_clip]]
=== Minimizing Fusion Transcripts Derived from Gene Dense Genomes (using --jaccard_clip)  ===

If your transcriptome RNA-seq data are derived from a gene-dense compact genome, such as from fungal genomes, where transcripts may often overlap in UTR regions, you can minimize fusion transcripts by leveraging the *--jaccard_clip* option if you have paired reads.  Trinity will examine the consistency of read pairings and fragment transcripts at positions that have little read-pairing support.  In expansive genomes of vertebrates and plants, this is unnecessary and not recommended.  In compact fungal genomes, it is highly recommended.  In addition to requiring paired reads, you must also have the http://bowtie-bio.sourceforge.net/index.shtml[Bowtie] short read aligner installed.  As part of this analysis, reads are aligned to the Inchworm contigs using Bowtie, and read pairings are examined across the Inchworm contigs, and contigs are clipped at positions of low pairing support.  These clipped Inchworm contigs are then fed into Chrysalis for downstream processing.  

Note, by using strand-specific RNA-Seq data alone, you should greatly mitigate the incorrect fusion of minimally overlapping transcripts.

[[genome_guided]]
=== Genome-guided Trinity  ===

If a genome sequence is available, Trinity offers a method whereby reads are first aligned to the genome, partitioned according to locus, followed by de novo transcriptome assembly at each locus.

Users must provide read alignments to Trinity as a coordinate-sorted bam file.  Use http://research-pub.gene.com/gmap/[gsnap], http://ccb.jhu.edu/software/tophat/index.shtml[tophat], https://github.com/alexdobin/STAR[STAR] or other favorite RNA-Seq read alignment tool to generate the bam file, and be sure it's coordinate sorted by running 'samtools sort' on it.

To run Genome-guided Trinity and have Trinity execute GSNAP to align the reads, run Trinity like so:

  Trinity --genome_guided_bam rnaseq.coordSorted.bam --genome_guided_max_intron 10000 --max_memory 10G --CPU 10 

Of course, use a maximum intron length that makes most sense given your targeted organism.

Be sure to include additional options such as '--SS_lib_type' and '--jaccard_clip' where appropriate.  If quality trimming or normalization are indicated, these processes will be performed prior to aligning the reads to the genome.

If you specify --grid_conf <string>, then the commands in this second phase will be executed in parallel on your compute farm, using LSF, SGE, or other supported method.  Otherwise, these commands will be executed locally using our Parafly parallel command processor, throttled at --CPU number of parallel processes.

[[genome_annotation]]
=== Comprehensive transcriptome-based genome annotation using Trinity and PASA ===

The Trinity-reconstructed transcripts can be used to annotate genomes using PASA.  Documentation for this is provided on the PASA website under http://pasapipeline.github.io/#A_ComprehensiveTranscriptome[Build a Comprehensive Transcriptome Database Using Genome-guided and De novo RNA-Seq Assembly] link.

[[trinity_output]]
== Output of Trinity ==

When Trinity completes, it will create a 'Trinity.fasta' output file in the 'trinity_out_dir/' output directory (or output directory you specify).  

Trinity groups transcripts into clusters based on shared sequence content. Such a transcript cluster is very loosely referred to as a 'gene'. This information is encoded in the Trinity fasta accession.  An example Fasta entry for one of the transcripts is formatted like so:

 >TR1000|c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]
 AATCTTTTTTGGTATTGGCAGTACTGTGCTCTGGGTAGTGATTAGGGCAAAAGAAGACAC
 ACAATAAAGAACCAGGTGTTAGACGTCAGCAAGTCAAGGCCTTGGTTCTCAGCAGACAGA
 AGACAGCCCTTCTCAATCCTCATCCCTTCCCTGAACAGACATGTCTTCTGCAAGCTTCTC
 CAAGTCAGTTGTTCACAGGAACATCATCAGAATAAATTTGAAATTATGATTAGTATCTGA
 TAAAGCA

The accession encodes the Trinity 'gene' and 'isoform' information. In the example above, the accession 'TR1000|c115_g5_i1' indicates Trinity read cluster 'TR1000|c115', gene 'g5', and isoform 'i1'.  Because a given run of trinity involves many many clusters of reads, each of which are assembled separately, and because the 'gene' numberings are unique within a given processed read cluster, the 'gene' identifier should be considered an aggregate of the read cluster and corresponding gene identifier, which in this case would be 'TR1000|c115_g5'.

So, in summary, the above example corresponds to 'gene id: TR1000|c115_g5' encoding 'isoform id: TR1000|c115_g5_i1'.


Obtain basic stats for the number of 'genes' and 'isoforms' and contiguity of the assembly by running:

  % $TRINITY_HOME/util/TrinityStats.pl trinity_out_dir/Trinity.fasta

with output (example from assembling our 10M Schizosaccharoymyces pombe data set):

 ################################
 ## Counts of transcripts, etc.
 ################################
 Total trinity 'genes':  8645
 Total trinity transcripts:  9398
 Percent GC: 37.59
 
 ########################################
 Stats based on ALL transcript contigs:
 ######################################## 

    Contig N10: 3838
    Contig N20: 3124
    Contig N30: 2629
    Contig N40: 2243
    Contig N50: 1936

    Median contig length: 984
    Average contig: 1251.23
    Total assembled bases: 11759032


 #####################################################
 ## Stats based on ONLY LONGEST ISOFORM per 'GENE':
 #####################################################

    Contig N10: 3848
    Contig N20: 3124
    Contig N30: 2630
    Contig N40: 2250
    Contig N50: 1937

    Median contig length: 942
    Average contig: 1227.97
    Total assembled bases: 10615785



[[compute_requirements]]
== Hardware and Configuration Requirements ==

The Inchworm and Chrysalis steps can be memory intensive.  A basic recommendation is to have ~1G of RAM per ~1M pairs of Illumina reads. Simpler transcriptomes (lower eukaryotes) require less memory than more complex transcriptomes such as from vertebrates.  

If you are able to run the entire Trinity process on a single high-memory multi-core server, indicate the number of butterfly processes to run in parallel by the --CPU parameter. 

Our experience is that the entire process can require ~1/2 hour to one hour per million pairs of reads in the current implementation (see link:trinity_faq.html[FAQ]).  We're striving to improve upon both memory and time requirements.


If you do not have direct access to a high memory machine (typically having 256G or 512G of RAM), consider <<RunElsewhere, running Trinity on one of the externally available resources>>.


[[monitoring_trinity]]
== Monitoring the Progress of Trinity ==
Since Trinity can easily take several days to complete, it is useful to be able to monitor the process and to know at which stage (Inchworm, Chrysalis, Butterfly) Trinity is currently at.  There are a few general ways to do this:

- by running 'top', you'll be able to see which Trinity process is running and how much memory is being consumed.
- other downstream process will generate standard output.  Be sure to capture 'stdout' and 'stderr' when you run the Trinity script.  The format for capturing both stdout and stderr depends on your SHELL.  Figure out what shell you have by running:

      env | grep SHELL

    Using tcsh:

         Trinity ... opts ... > & run.log &

    Using bash:

        Trinity ... opts ... > run.log 2>&1 &

Note, under bash, to prevent the background process from being terminated once you close the shell, type 'exit' to leave the shell, or explore alternatives such as http://www.serverwatch.com/tutorials/article.php/3935306/Detach-Processes-With-Disown-and-Nohup.htm[nohup, disown, or screen].

You can then 'tail -f run.log' to follow the progress of the Trinity throughout the various stages.


[[sample_data]]
== Running Trinity on Sample Data ==

The Trinity software distribution includes sample data in the 'sample_data/test_Trinity_Assembly/' directory. Simply run the included 'runMe.sh' shell script to execute the Trinity assembly process with provided paired strand-specific Illumina data derived from mouse.  Running Trinity on the sample data requires <~2G of RAM and should run on an ordinary desktop/laptop computer.  Run as 'runMe.sh 1' to execute downstream analysis steps, including bowtie read alignment and RSEM-based abundance estimation, as described below.


[[Downstream_analyses]]
== Downstream Analyses ==

The following downstream analyses are supported as part of Trinity:

- link:analysis/abundance_estimation.html[Abundance estimation using RSEM or eXpress, and visualization using IGV].
- link:analysis/diff_expression_analysis.html[Using EdgeR and Bioconductor for analyzing differentially expressed transcripts].
- http://transdecoder.github.io[Extract likely protein-coding regions from Trinity transcripts using TransDecoder].
- http://trinotate.github.io[Functionally annotate transcripts and coding regions with Trinotate].
- link:analysis/full_length_transcript_analysis.html[Full-length transcript analysis for model and non-model transcriptomes]

[[advanced_guide]]
== Want to know more? ==

Visit the link:advanced_trinity_guide.html[Advanced Guide to Trinity] for more information regarding Trinity behavior, intermediate data files, and file formats.

[[faq]]
== Frequently Asked Questions ==

Visit the link:trinity_faq.html[Trinity FAQ] page.

[[trinity_tidbits]]
== Trinity Tidbits ==

- Trinity made the cover of the http://www.nature.com/nbt/journal/v29/n7/index.html[July 2011 NBT issue]. The Broad Institute's http://www.broadinstitute.org/blog/suite-tools-takes-flight[blog] has a story on how the Trinity project came together. Nir Friedman, one of the project PIs, has a http://nirfriedmanlab.blogspot.com/2011/07/behind-cover.html[blog entry] describing the developmental process underlying the NBT cover design.

- Trinity was shown to be the leading de novo transcriptome assembly tool as part of the http://www.the-dream-project.org/challanges/dream6-alternative-splicing-challenge[DREAM6 Alt-Splicing Challenge 2011]. Results were posted http://www.the-dream-project.org/result/alternative-splicing[here].  

- http://scholar.google.com/scholar?oi=bibs&hl=en&cites=14735674943942667509[Google Scholar] shows how Trinity is being used by the community.

[[contact_us]]
== Contact Us ==

Questions, suggestions, comments, etc?

Join and add discussions at the Trinityrnaseq-users Google group: https://groups.google.com/forum/\#!forum/trinityrnaseq-users[https://groups.google.com/forum/#!forum/trinityrnaseq-users].


[[referencing_trinity]]
== Referencing Trinity ==

Trinity can be referenced as:

- Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N,
di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A.
Full-length transcriptome assembly from RNA-seq data without a reference genome. 
http://www.nature.com/nbt/journal/vaop/ncurrent/abs/nbt.1883.html[Nat Biotechnol. 2011 May 15;29(7):644-52]. doi: 10.1038/nbt.1883. 
http://www.ncbi.nlm.nih.gov/pubmed/21572440[PubMed PMID: 21572440].

Protocol for using Trinity for de novo transcriptome assembly and downstream analyses:

- Haas BJ, Papanicolaou A, Yassour M, Grabherr M, Blood PD, Bowden J, Couger MB,
Eccles D, Li B, Lieber M, Macmanes MD, Ott M, Orvis J, Pochet N, Strozzi F, Weeks
N, Westerman R, William T, Dewey CN, Henschel R, Leduc RD, Friedman N, Regev A.
De novo transcript sequence reconstruction from RNA-seq using the Trinity
platform for reference generation and analysis. http://www.nature.com/nprot/journal/v8/n8/full/nprot.2013.084.html[Nat Protoc. 2013 Aug;8(8):1494-512.] http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/[Open Access in PMC] doi: 10.1038/nprot.2013.084. Epub 2013 Jul 11. PubMed PMID:
23845962.


Performance tuning of Trinity is described in:

- Henschel R, Lieber M, Wu L, Nista, PM, Haas BJ, LeDuc R.  Trinity RNA-Seq assembler performance optimization. XSEDE 2012 Proceedings of the 1st Conference of the Extreme Science and Engineering Discovery Environment: Bridging from the eXtreme to the campus and beyond. http://dx.doi.org/10.1145/2335755.2335842[ISBN: 978-1-4503-1602-6 doi>10.1145/2335755.2335842].

A full list of references including Trinity, RSEM, and additional tools leveraged by Trinity can be obtained by running 'Trinity --cite'.
