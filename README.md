# Phastats

A gene sequence statistics tool

In comparison to Fastqc, our goal is to make a tool that takes in an argument fastq file and outputs a html file containing similar information. In running our Sequence Statistics tool, the metrics of the resulting html file would include quality distribution, length distribution, and general, useful counts including mean read length and mean quantity of the reads.

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
