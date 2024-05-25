# Phastats

A gene sequence statistics tool

In comparison to Fastqc, our goal is to make a tool that takes in an argument fastq file and outputs a html file containing similar information. In running our Sequence Statistics tool, the metrics of the resulting html file would include quality distribution, length distribution, and general, useful counts including mean read length and mean quantity of the reads.

### Emily:

My file calculates the per-base sequence content at each position in the read. For simplicity, it does not take in any command line arguments and the fasta file that I was working with is hard coded into my py file. It can be run by typing ```python Emily.py```` and a graph and table should display.
