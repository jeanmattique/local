#Running Trinity

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



Trinity performs best with strand-specific data, in which case sense and antisense transcripts can be resolved.  For protocols on strand-specific RNA-Seq, see: http://www.ncbi.nlm.nih.gov/pubmed/21943893[Borodina T, Adjaye J, Sultan M. A strand-specific library preparation protocol for RNA sequencing. Methods Enzymol. 2011;500:79-98. PubMed PMID: 21943893].


If you have strand-specific data, specify the library type.  There are four library types:

- Paired reads:
    * *RF*: first read (/1) of fragment pair is sequenced as anti-sense (reverse(*R*)), and second read (/2) is in the sense strand (forward(*F*)); typical of the dUTP/UDG sequencing method.
    * *FR*: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the antisense strand (reverse)

- Unpaired (single) reads:
    * *F*: the single read is in the sense (forward) orientation
    * *R*: the single read is in the antisense (reverse) orientation

By setting the *--SS_lib_type* parameter to one of the above, you are indicating that the reads are strand-specific.  By default, reads are treated as not strand-specific.

![strand specific specification](https://raw.githubusercontent.com/wiki/trinityrnaseq/trinityrnaseq/images/strand_specificity.jpg)

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