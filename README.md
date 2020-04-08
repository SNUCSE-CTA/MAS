# Fast String Matching for DNA Sequences  [Theoretical Computer Science 2019]
Cheol Ryu, Thierry Lecroq, and Kunsoo Park

## MAS, TMAS, QMAS
Maximal Average Shift (MAS) algorithm

Tuned Maximal Average Shift (TMAS) algorithm

q-gram Maximal Average Shift (QMAS) algorithm

## Introduction
In this paper we propose the Maximal Average Shift (MAS) algorithm that finds a pattern scan order that maximizes the average shift length. We also present two extensions of MAS: one improves the scan speed of MAS by using the scan result of the previous window, and the other improves the running time of MAS by using $q$-grams. These algorithms show better average performances in scan speed than previous string matching algorithms for DNA sequences.

## Source Code Information
This source code is an implementation of \[1\], fast string matching algorithms for DNA Sequences.

The source code is implemented based on SMART \[2\]. 

(SMART is available at https://github.com/smart-tool/smart)

## Datasets
Human chromosomes 20 downloaded from the 1000 Genomes Project website \[3\].

## Run

$ ./mas (pattern) (pattern length) (text) (text length)

$ ./mas GCAGAGAG 8 GCATCGCAGAGAGTATACAGTACG 24

(∑ = A,C,G,T)

## References
[1] C. Ryu, T. Lecroq, and K. Park, Fast String Matching for DNA Sequences, Theoretical Computer Science 812 (2020) 137-148.

[2] S. Faro, T. Lecroq, S. Borzi, S. D. Mauro, A. Maggio, The string matching algorithms research tool, in: Proceedings of the Prague Stringology Conference 2016, Prague Stringology Club, 2016, pp. 99-113.

[3] The 1000 Genomes Project Consortium, The 1000 Genomes Project, http://www.1000genomes.org, 2010.


