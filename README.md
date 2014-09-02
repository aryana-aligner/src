src
===

The source code of the Aryana next generation sequencing (NGS) aligner.

Installation
============

Download Aryana from github and run "make"

Running Aryana
==============

Alignment of single end reads

<code>
./aryana [-x human_g1k_v37.fasta] [-U reads.fastq]

	-p INT	Number of threads

	-P INT	Number of potential candidates to select for more precise alignment

	--seed INT	Fixed length for the seeds
</code>

Alignment of paired end reads

<code>
./aryana [-x human_g1k_v37.fasta] [-1 pair1.fastq] [-2 pair2.fastq] --{fr, ff, rf} -I min -X max

	--fr/--ff/--rf	Refers to orientation of the pairs. /forward-reverse/forward-forward/reverse-forward

	-I INT	Minimum insert size

	-X INT	Maximum insert size
</code>

