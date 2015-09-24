#Output of Trinity Assembly

When Trinity completes, it will create a 'Trinity.fasta' output file in the 'trinity_out_dir/' output directory (or output directory you specify).  

Trinity groups transcripts into clusters based on shared sequence content. Such a transcript cluster is very loosely referred to as a 'gene'. This information is encoded in the Trinity fasta accession.  An example Fasta entry for one of the transcripts is formatted like so:

     >TRINITY_DN1000|c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]
     AATCTTTTTTGGTATTGGCAGTACTGTGCTCTGGGTAGTGATTAGGGCAAAAGAAGACAC
     ACAATAAAGAACCAGGTGTTAGACGTCAGCAAGTCAAGGCCTTGGTTCTCAGCAGACAGA
     AGACAGCCCTTCTCAATCCTCATCCCTTCCCTGAACAGACATGTCTTCTGCAAGCTTCTC
     CAAGTCAGTTGTTCACAGGAACATCATCAGAATAAATTTGAAATTATGATTAGTATCTGA
     TAAAGCA

The accession encodes the Trinity 'gene' and 'isoform' information. In the example above, the accession 'TRINITY_DN1000|c115_g5_i1' indicates Trinity read cluster 'TRINITY_DN1000|c115', gene 'g5', and isoform 'i1'.  Because a given run of trinity involves many many clusters of reads, each of which are assembled separately, and because the 'gene' numberings are unique within a given processed read cluster, the 'gene' identifier should be considered an aggregate of the read cluster and corresponding gene identifier, which in this case would be 'TRINITY_DN1000|c115_g5'.

So, in summary, the above example corresponds to 'gene id: TRINITY_DN1000|c115_g5' encoding 'isoform id: TRINITY_DN1000|c115_g5_i1'.


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
