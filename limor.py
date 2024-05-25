#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt

def read_lengths_and_quality_scores_from_fastq(file_path):
    lengths = []
    quality_scores = []
    with open(file_path, "r") as handle:
        while True:
            header = handle.readline().strip()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline().strip()
            quality = handle.readline().strip()
            lengths.append(len(sequence))
            quality_scores.extend([ord(q) - 33 for q in quality])  # Convert ASCII quality scores to Phred scores
    return lengths, quality_scores

def plot_length_distribution(lengths, output_file):
    plt.figure(figsize=(9, 5))
    plt.hist(lengths, bins=50, color='#0059b3', alpha=0.5, edgecolor='black')  # Blue color
    plt.title('Length Distribution of Sequences', fontsize=16, fontweight='bold')
    plt.xlabel('Sequence Length', fontsize=13)
    plt.ylabel('Frequency', fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='dotted')
    plt.savefig(output_file + "_length_distribution.png", dpi=300)  # Higher resolution
    plt.show()

def plot_quality_distribution(quality_scores, output_file):
    plt.figure(figsize=(9, 5))
    plt.hist(quality_scores, bins=range(0, 41), color='#008000', alpha=0.5, edgecolor='black')  # Green color
    plt.title('Quality Distribution of Sequences', fontsize=16, fontweight='bold')
    plt.xlabel('Quality Score (Phred Score)', fontsize=13)
    plt.ylabel('Frequency', fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='dotted')
    plt.savefig(output_file + "_quality_distribution.png", dpi=300)  # Higher resolution
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot length and quality distribution of sequences in a FASTQ file.')
    parser.add_argument('input_file', type=str, help='Path to the input FASTQ file.')
    parser.add_argument('output_file', type=str, help='Prefix for the output plot filenames.')
    args = parser.parse_args()

    lengths, quality_scores = read_lengths_and_quality_scores_from_fastq(args.input_file)
    plot_length_distribution(lengths, args.output_file)
    plot_quality_distribution(quality_scores, args.output_file)

if __name__ == "__main__":
    main()
