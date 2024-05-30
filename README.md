# Phastats

A gene sequence statistics tool

In comparison to Fastqc, our goal is to make a tool that takes in an argument fastq file and outputs a html file containing similar information. In running our Sequence Statistics tool, the metrics of the resulting html file would include N50 score, quality distribution, length distribution, and general, useful counts including mean read length and mean quantity of the reads.

## Prerequisites

`phastats` requires:

- Python 3.xx
- Python packages:
  - pandas
  - numpy
  - scipy
  - matplotlib
  - argparse

If you do not have these packages, you can intall them by running this command:

```
pip install pandas numpy scipy matplotlib
```

If you do not have root access, you can run the command above with the additional `--user` option to install locally:

```
pip install --user pandas numpy scipy matplotlib
```

## Implementation

This Python script conducts comprehensive analysis on FASTQ datasets and generates informative plots. Here's a detailed summary of how the code is implemented:

First, Phastats imports its necessary libraries as shown above for functions including parsing, numerical operations and statistical functions, data manipulation, and other operations.

Using a .fq dataset as an argument, Phastats parses command-line arguments to get FASTQ file information such as sequence lengths, quality scores, and GC content; Phastats reads four lines at a time for each record. The following functions 
```getandPlotLengths``` 
reads the sequences, ```getandPlotQuality``` quality lines in each record, and parse_fastq do this to extract this information.

In addition to calculating the number of sequences in addition to the length of each sequence, Phastats converts the ASCII quality scores to Phred scores and 'G' and 'C' nucleotides for analysis regarding quality distribution and GC content.

## Running phastats




## Benchmarking

for when we finish benchmarking it lol


There are various reasons why the outputs of our phastats tool contrasts from the FastQC tool. To start, our Phastats tool computes the N50 of the .fq dataset, while FastQC does not. This is useful, for the N50 score provides another measure of the quality of a genome assembly in addition to the information displayed in the resulting html file.
Another reason why these tools differ could be that FastQC was coded in Java, while Phastats was created in Python. Because of this, Phastats contains simple, readable syntax, allowing developers to express concepts in fewer lines of code. 
Finally, we have implemented a more modern design for our Phastats HTML file, which may enhance readability compared to FastQC.



