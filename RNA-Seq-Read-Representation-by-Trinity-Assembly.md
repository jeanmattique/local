#Assessing the Read Content of the Transcriptome Assembly

Assembled transcripts might not always fully represent properly paired-end reads, as some transcripts may be fragmented or short and only one fragment read of a pair may align.  Simply aligning reads to your transcriptome assembly using bowtie or STAR will only capture the properly paired reads.  To assess the read composition of our assembly, we want to capture and count all reads that map to our assembled transcripts, including the properly paired and those that are not.

In order to comprehensively capture read alignments, we run the process below.  Bowtie is used to align each fragment end to the transcriptome assembly separately. Subsequently, the read pairs are grouped into properly paired reads where possible, and those reads that do not map as properly paired are still retained.  

      $TRINITY_HOME/util/bowtie_PE_separate_then_join.pl --seqType fq \
                  --left left.fq --right right.fq \
                  --target Trinity.fasta --aligner bowtie \
                  -- -p 4 --all --best --strata -m 300  # following -- are params that get tacked onto the bowtie command.


As usual, if you have strand-specific RNA-Seq data, indicate this with the '--SS_lib_type' parameter, and put this parameter before the '--' above, since all the parameters after '--' are applied to the bowtie aligner.

An output directory 'bowtie_out' is created and should include the files:

     bowtie_out.nameSorted.bam  : alignments sorted by read name
     bowtie_out.coordSorted.bam : alignments sorted by coordinate.


To get alignment statistics, run the following on the name-sorted bam file:

       $TRINITY_HOME/util/SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.bam

     #read_type  count   pct
     proper_pairs    47042   83.59  (left and right reads align to the same transcript)
     improper_pairs  6824    12.13  (left and right reads align, but to different transcripts)
     left_only   1300    2.31 
     right_only  1110    1.97
 
     Total aligned reads: 56276  (counting individual reads of pairs, each read counts only once).



A typical Trinity transcriptome assembly will have the vast majority of all reads mapping back to the assembly, and ~70-80% of the mapped reads found mapped as proper pairs.

