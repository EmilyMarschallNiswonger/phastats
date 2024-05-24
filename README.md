# Phastats

A gene sequence statistics tool

In comparison to Fastqc, our goal is to make a tool that takes in an argument fastq file and outputs a html file containing similar information. In running our Sequence Statistics tool, the metrics of the resulting html file would include quality distribution, length distribution, and general, useful counts including mean read length and mean quantity of the reads.
Somto: For portion of the progress report I was able to create a python script that first gives a table of the Basic Statisics of a fastqc files (Filesname, File type, Encoding, Total Sequences, Sequences flagged as poor quality, Sequence Length, and %GC conctent). The other portion give you a bar/line graph that highlights the per sequence GC content. All These are runnable using the code "python Somto.py samples/data1.fq
