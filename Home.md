#RNA-Seq De novo Assembly Using Trinity

![TrinityCompositeLogo](https://raw.githubusercontent.com/wiki/trinityrnaseq/trinityrnaseq/images/TrinityCompositeLogo.png)

##Quick Guide for the Impatient

Trinity assembles transcript sequences from Illumina RNA-Seq data.

Download Trinity [here](https://github.com/trinityrnaseq/trinityrnaseq/releases).

Build Trinity by typing 'make' in the base installation directory.

Assemble RNA-Seq data like so:

     Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G 

Find assembled transcripts as:  'trinity_out_dir/Trinity.fasta'

##Intro to Trinity

Trinity, developed at the [Broad Institute](http://www.broadinstitute.org) and the [Hebrew University of Jerusalem] (http://www.cs.huji.ac.il), represents a novel method for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data. Trinity combines three independent software modules: Inchworm, Chrysalis, and Butterfly, applied sequentially to process large volumes of RNA-seq reads. Trinity partitions the sequence data into many individual de Bruijn graphs, each representing the transcriptional complexity at at a given gene or locus, and then processes each graph independently to extract full-length splicing isoforms and to tease apart transcripts derived from paralogous genes.  Briefly, the process works like so:

- *Inchworm* assembles the RNA-seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.

- *Chrysalis* clusters the Inchworm contigs into clusters and constructs complete de Bruijn graphs for each cluster.  Each cluster represents the full transcriptonal complexity for a given gene (or sets of genes that share sequences in common).  Chrysalis then partitions the full read set among these disjoint graphs.

- *Butterfly* then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.

Trinity was published in [Nature Biotechnology](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571712/).  Our protocol for transcriptome assembly and downstream analysis is published in [Nature Protocols] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/), although we always have the most current instructional material available here at the Trinity website.

The Trinity software package can be downloaded [here on GitHub](https://github.com/trinityrnaseq/trinityrnaseq/releases). Legacy versions (pre-2015) are still available at [our Sourceforge Trinity software archive](http://sourceforge.net/projects/trinityrnaseq/files/PREV_CONTENTS/previous_releases/).

[Runtime and transcript reconstruction performance stats](http://trinityrnaseq.github.io/performance/) are available for current and previous releases.

[Screenast videos](http://www.broadinstitute.org/partnerships/education/broade/trinity-screencast) are available to introduce you to Trinity and its various components. Also, hands-on tutorials for Trinity and Tuxedo are available as part of our link:workshop/rnaseq_workshop.html[RNA-Seq Workshop].




== Table of Contents ==


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
