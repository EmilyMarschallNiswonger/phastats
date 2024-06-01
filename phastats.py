#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import os
import pandas as pd

def getandPlotLengths(sequenceList):
    lengths = [len(sequence) for sequence in sequenceList]
    return plot_length_distribution(lengths)

def getandPlotQuality(qualityList):
    quality_scores = [ord(q) - 33 for q in qualityList[0]]  # Convert ASCII quality scores to Phred scores
    return plot_quality_distribution(quality_scores)


def plot_length_distribution(lengths):
    plt.figure(figsize=(9, 5))
    plt.hist(lengths, bins=50, color='#0059b3', alpha=0.5, edgecolor='black')  # Blue color
    plt.title('Length Distribution of Sequences', fontsize=16, fontweight='bold')
    plt.xlabel('Sequence Length', fontsize=13)
    plt.ylabel('Frequency', fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='dotted')
    plt.savefig("plots/length_distribution.png", dpi=300)  # Higher resolution

def plot_quality_distribution(quality_scores):
    plt.figure(figsize=(9, 5))
    plt.hist(quality_scores, bins=range(0, 41), color='#008000', alpha=0.5, edgecolor='black')  # Green color
    plt.title('Quality Distribution of Bases', fontsize=16, fontweight='bold')
    plt.xlabel('Quality Score (Phred Score)', fontsize=13)
    plt.ylabel('Frequency', fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='dotted')
    plt.savefig("plots/quality-distribution.png", dpi=300)

def getandPlotGCDistribution(gc_count, total_length, gc_contents):
    gc_content = calculate_gc_content(gc_count, total_length)

    # print_statistics(file_path, total_sequences, poor_quality_sequences, total_length, gc_content)
    plot_gc_distribution(gc_contents)

def parse_fastq(file):
    total_sequences = 0
    poor_quality_sequences = 0
    total_length = 0
    gc_count = 0
    gc_contents = []
    sequences = []
    qualitylines = []
    i = 0;

    with open(file, 'r') as f:
        while True:
            # Read the four lines of a FASTQ record, create a list of sequences and quality lines
            identifier = f.readline().strip()
            if not identifier:
                break
            sequences.append(f.readline().strip())
            plus = f.readline().strip()
            qualitylines.append(f.readline().strip())

            total_sequences += 1
            total_length += len(sequences[i])
            gc_count += sequences[i].count('G') + sequences[i].count('C')
            
            # Calculate GC content for the sequence
            sequence_gc_count = sequences[i].count('G') + sequences[i].count('C')
            sequence_gc_content = (sequence_gc_count / len(sequences[i])) * 100
            gc_contents.append(sequence_gc_content)

            # Check if the sequence is flagged as poor quality
            if any(char < '!' for char in qualitylines[i]):  # Placeholder for poor quality check
                poor_quality_sequences += 1
            i += 1

    return total_sequences, poor_quality_sequences, total_length, sequences, qualitylines, gc_count, gc_contents

def calculate_gc_content(gc_count, total_length):
    return (gc_count / total_length) * 100 if total_length > 0 else 0

def print_statistics(filename, total_sequences, poor_quality_sequences, total_length, gc_content):
    average_length = total_length / total_sequences if total_sequences > 0 else 0

    print(f"{'Measure':<30}{'Value':<15}")
    print("-" * 45)
    print(f"{'Filename':<30}{filename:<15}")
    print(f"{'File type':<30}{'FASTQ':<15}")
    print(f"{'Encoding':<30}{'Sanger / Illumina 1.8+':<15}")  # Assumed encoding
    print(f"{'Total sequences':<30}{total_sequences:<15}")
    print(f"{'Sequences flagged as poor quality':<30}{poor_quality_sequences:<15}")
    print(f"{'Average sequence length':<30}{average_length:<15.2f}")
    print(f"{'GC content (%)':<30}{gc_content:<15.2f}")

def plot_gc_distribution(gc_contents):
    plt.figure(figsize=(10, 6))

    # Histogram for GC content per read
    count, bins, ignored = plt.hist(gc_contents, bins=50, density=True, alpha=0.6, color='r', edgecolor='black', label='GC count per read')

    # Fit a normal distribution to the data
    mu, std = np.mean(gc_contents), np.std(gc_contents)
    p = stats.norm.pdf(bins, mu, std)

    # Plot the theoretical normal distribution
    plt.plot(bins, p, 'b-', linewidth=2, label='Theoretical Distribution')

    plt.xlabel('Mean GC content (%)')
    plt.ylabel('Frequency')
    plt.title('GC distribution over all sequences')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig("plots/gc-distribution.png", dpi=300)
    plt.clf()

def getandPlotPerBaseSequenceContent(lines):
    # Read the FASTQ file and calculate the percent of each base per read
    df = pd.DataFrame(columns=['%A', '%G', '%C', '%T'])
    max_length = max(len(line) for line in lines)
    for i in range(max_length):
        A = sum(1 for line in lines if i < len(line) and line[i] == 'A')
        G = sum(1 for line in lines if i < len(line) and line[i] == 'G')
        C = sum(1 for line in lines if i < len(line) and line[i] == 'C')
        T = sum(1 for line in lines if i < len(line) and line[i] == 'T')
        total = A + G + C + T
        if total > 0:  # avoid division by zero
            df.loc[i] = [A/total, G/total, C/total, T/total]
    # plot the percent of each base per read. y axis is percent, and 
    # x axis is read position, or the row numnber
    # write the plot to an html file
    plt.plot(df['%A'], label='%A')  
    plt.plot(df['%G'], label='%G')
    plt.plot(df['%C'], label='%C')
    plt.plot(df['%T'], label='%T')
    plt.xlabel('Read Position')
    plt.ylabel('Percent')
    plt.title('Per-Base Sequence Content')
    plt.legend()
    plt.savefig("plots/Per-base_sequnce_content.png", dpi=300)
    plt.clf()

def compute_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    cumulative_length = sum(sorted_lengths)
    half_length = cumulative_length / 2
    running_length = 0
    for length in sorted_lengths:
        running_length += length
        if running_length >= half_length:
            return length


def main():
    parser = argparse.ArgumentParser(description='Conduct fastq analysis and plot the results.')
    parser.add_argument('input_file', type=str, help='Path to the input FASTQ file.')
    # parser.add_argument('output_file', type=str, help='Prefix for the output plot filenames.')
    args = parser.parse_args()


    # Parse the FASTQ file and calculate statistics
    total_sequences, poor_quality_sequences, total_length, sequences, qualitylines, gc_count, gc_contents = parse_fastq(args.input_file)
    
    # Compute N50
    n50_value = compute_n50(lengths)

    # Generage plot for length and quality distribution
    getandPlotLengths(sequences)
    getandPlotQuality(qualitylines)
    
    # Generate plot for gc distribution. Not printing statistics for now
    getandPlotGCDistribution(gc_count, total_length, gc_contents)

    # Generate plot for perbase sequence content
    getandPlotPerBaseSequenceContent(sequences)




if __name__ == "__main__":
    main()