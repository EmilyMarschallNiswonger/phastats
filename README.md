# Phastats

A gene sequence statistics tool

In comparison to Fastqc, our goal is to make a tool that takes in an argument fastq file and outputs a html file containing similar information. In running our Sequence Statistics tool, the metrics of the resulting html file would include quality distribution, length distribution, and general, useful counts including mean read length and mean quantity of the reads.


### Emily:

My file calculates the per-base sequence content at each position in the read. For simplicity, it does not take in any command line arguments and the fasta file that I was working with is hard coded into my py file. It can be run by typing ```python Emily.py```` and a graph and table should display. I did not yet make the file executable

### Somto:
For portion of the progress report I was able to create a python script that first gives a table of the Basic Statisics of a fastqc files (Filesname, File type, Encoding, Total Sequences, Sequences flagged as poor quality, Sequence Length, and %GC conctent). The other portion give you a bar/line graph that highlights the per sequence GC content. All These are runnable using the code "python Somto.py samples/data1.fq

### Li:

My file “limor.py” calculates the per sequence quality and length distribution of an argument FASTQ file. This file is an executable; this is designed to be run from the command line, so it could be run using command, “./limor.py samples/data1.fq output_plot.png”. This example command takes in sample file data1.fq and outputs two pngs, output_plot_quality_distribution.png and output_plot_length_distribution.png.