






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

## Running phastats




## Implementation

This Python script conducts comprehensive analysis on FASTQ datasets and generates informative plots. Here's a detailed summary of how the code is implemented:

First, Phastats imports its necessary libraries as shown above for functions including parsing, numerical operations and statistical functions, data manipulation, and other operations.

Using a .fq dataset as an argument, Phastats parses command-line arguments to get FASTQ file information such as sequence lengths, quality scores, and GC content; Phastats reads four lines at a time per sequence. 


`getLengthAndQuality` reads a FASTQ file and extracts sequence lengths and quality scores. Quality scores are converted from ASCII characters to Phred scores.

`plot_length_distribution` creates and saves a histogram of sequence lengths to the final html file.

`plot_quality_distribution` creates and saves a histogram of quality scores to the final html file.

`getLengthQualityDistribution` calls `getLengthAndQuality` and the respective plotting functions `plot_length_distribution` and `plot_quality_distribution` to generate the length and quality distribution plots.

`parse_fastq` parses the FASTQ file to extract the total number of sequences and poor quality sequences, total length, GC count, and GC content for each sequence.

`calculate_gc_content` calculates the overall GC content percentage by dividing the GC count by the total length.

`plot_gc_distribution` creates and saves a histogram of GC content per sequence to the final html file.

`getGCDistribution` calls `parse_fastq` and `plot_gc_distribution` to generate the GC content plot.

`getPerBaseSequenceContent` calculates and plots the percentage of each base (A, G, C, T) at each position in each sequence.

`compute_n50` computes the N50 value, a measure of the quality of genome assemblies, by sorting sequence lengths and finding the length at which 50% of the total sequence length is contained.

`print_statistics` outputs various metrics for the provided FASTQ file such as filename, file type encoding, total count of sequences, sequences flagged as poor quality, average sequence length and percentage GC content. 




## Benchmarking

for when we finish benchmarking it lol


There are various reasons why the outputs of our phastats tool contrasts from the FastQC tool. To start, our Phastats tool computes the N50 of the .fq dataset, while FastQC does not. This is useful, for the N50 score provides another measure of the quality of a genome assembly in addition to the information displayed in the resulting html file.
Another reason why these tools differ could be that FastQC was coded in Java, while Phastats was created in Python. Because of this, Phastats contains simple, readable syntax, allowing developers to express concepts in fewer lines of code. 
Finally, we have implemented a more modern design for our Phastats HTML file, which may enhance readability compared to FastQC.





