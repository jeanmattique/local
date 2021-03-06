# Transcriptome Assembly Quality Assessment

Once your assembly is complete, you'll want to know how 'good' it is, and you might want to compare the quality of the assembly to similar assemblies generated by alternative assemblers, or having run an assembly with different parameters.

There are some general ways to characterize the quality of your assembly:

*  Examine the [RNA-Seq read representation of the assembly](RNA-Seq-Read-Representation-by-Trinity-Assembly). Ideally, at least ~80% of your input RNA-Seq reads are represented by your transcriptome assembly.  The remaining unassembled reads likely corresponds to lowly expressed transcripts with insufficient coverage to enable assembly, or are low quality or aberrant reads.

*   Examine the [representation of full-length reconstructed protein-coding genes](Counting-Full-Length-Trinity-Transcripts), by searching the assembled transcripts against a database of known protein sequences.

*   Compute the [E90N50 transcript contig length](Transcriptome-Contig-Nx-and-ExN50-stats) - the contig N50 value based on the set of transcripts representing 90% of the expression data.

*   Compute [DETONATE](http://deweylab.biostat.wisc.edu/detonate/) scores.  DETONATE provides a rigorous computational assessment of the quality of a transcriptome assembly.