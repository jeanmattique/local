#Full-length transcript analysis for model and non-model organisms using BLAST+

One metric for evaluating the quality of a transcriptome assembly is to examine the number of transcripts that were assembled that appear to be full-length or nearly full-length.  Such an analysis with a reference transcript set, such as from human or mouse, is relatively straightforward, since one can align the assembled transcripts to the reference transcripts and examine the length coverage.  For non-model organisms, no such reference transcript set is available. If a high quality annotation exists for a closely related organism, then one might compare the assembled transcripts to that closely related transcriptome to examine full-length coverage. In other cases, a more general analysis to perform is to align the assembled transcripts against all known proteins and to determine the number of unique top matching proteins that align across more than X% of its length.

Trinity supports these analyses using [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).  

Useful protein databases to search include [SwissProt](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) (<ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz>) and [TrEMBL](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz) (<ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz>).

To examine the extent of top-matching BLAST alignments, first run BLAST, and then run the included analysis script below:

     Build a blastable database:
     
     %  makeblastdb -in uniprot_sprot.fasta -dbtype prot
     
     
     Perform the blast search, reporting only the top alignment:
     
      % blastx -query Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 \
            -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6
