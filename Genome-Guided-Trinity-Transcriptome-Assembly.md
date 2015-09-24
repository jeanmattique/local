#Genome-guided Trinity De novo Transcriptome Assembly 

If a genome sequence is available, Trinity offers a method whereby reads are first aligned to the genome, partitioned according to locus, followed by de novo transcriptome assembly at each locus.

Users must provide read alignments to Trinity as a coordinate-sorted bam file.  Use [GSNAP](http://research-pub.gene.com/gmap/), [TopHat](http://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://github.com/alexdobin/STAR) or other favorite RNA-Seq read alignment tool to generate the bam file, and be sure it's coordinate sorted by running 'samtools sort' on it.

To run Genome-guided Trinity and have Trinity execute GSNAP to align the reads, run Trinity like so:

     Trinity --genome_guided_bam rnaseq.coordSorted.bam --genome_guided_max_intron 10000 --max_memory 10G --CPU 10 

Of course, use a maximum intron length that makes most sense given your targeted organism.

Be sure to include additional options such as '--SS_lib_type' and '--jaccard_clip' where appropriate.  If quality trimming  or normalization are desired, these processes should be performed prior to aligning the reads to the genome, as Trinity will only use the reads as they exist in the coordinate-sorted bam file provided to it.

If you specify --grid_conf <string>, then the commands in this second phase will be executed in parallel on your compute farm, using LSF, SGE, or other supported method.  Otherwise, these commands will be executed locally using our Parafly parallel command processor, throttled at --CPU number of parallel processes.