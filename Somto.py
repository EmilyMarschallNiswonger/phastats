#!/usr/bin/env python

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def parse_fastq(file):
    total_sequences = 0
    poor_quality_sequences = 0
    total_length = 0
    gc_count = 0
    gc_contents = []

    with open(file, 'r') as f:
        while True:
            # Read the four lines of a FASTQ record
            identifier = f.readline().strip()
            if not identifier:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            total_sequences += 1
            total_length += len(sequence)
            gc_count += sequence.count('G') + sequence.count('C')
            
            # Calculate GC content for the sequence
            sequence_gc_count = sequence.count('G') + sequence.count('C')
            sequence_gc_content = (sequence_gc_count / len(sequence)) * 100
            gc_contents.append(sequence_gc_content)

            # Check if the sequence is flagged as poor quality
            if any(char < '!' for char in quality):  # Placeholder for poor quality check
                poor_quality_sequences += 1

    return total_sequences, poor_quality_sequences, total_length, gc_count, gc_contents

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
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Calculate basic statistics for a FASTQ file.")
    parser.add_argument('input_file', type=str, help='Input FASTQ file')

    args = parser.parse_args()

    filename = os.path.basename(args.input_file)
    total_sequences, poor_quality_sequences, total_length, gc_count, gc_contents = parse_fastq(args.input_file)
    gc_content = calculate_gc_content(gc_count, total_length)

    print_statistics(filename, total_sequences, poor_quality_sequences, total_length, gc_content)
    plot_gc_distribution(gc_contents)

if __name__ == "__main__":
    main()