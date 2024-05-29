#!/usr/bin/env python 


"""
This file calculates the per-base sequence content at each position 
in the read
Takes in a fq file, outputs a graph and dataframe showing percent A
G C and T at each positon in the read
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import html
from pathlib import Path
filename = "samples/child1_big.fq"

def main():
    # read in file and create a list of strings, each string being every
    # 4th line in the file, starting from the second line
    with open(filename, 'r') as file:
        lines = file.readlines()
        lines = lines[1::4]
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
    print("Number of rows: ", df.shape[0])

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
    with (open("sequence.html", "w")) as f:
        f.write("<html><body><img src='plots/Per-base_sequnce_content.png' alt='plot'></body></html>")
    




if __name__ == "__main__":
    main()