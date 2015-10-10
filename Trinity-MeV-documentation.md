# Using MeV to navigate your Trinity DE genes/transcripts

Welcome to our introductory guide to using MeV to navigating your differentially expressed gene or transcript data!

First, be sure to have the [TM4 MeV](http://www.tm4.org/mev.html) software installed.

Next, follow the instructions below to upload your expression matrix into MeV.  We'll assume that you're starting with results generated from having run through the [DE analysis] step, and have a file such as 'diffExpr.P0.001_C2.matrix.log2.dat' in your workspace, as generated from 'analyze_diff_expr.pl' in the DE analysis step.  This matrix contains the expression values for those features that are at least 2^2 (or 4) fold DE with an FDR <= 0.001.  Note, as you can see from the filename, this matrix has already been log2 transformed.
