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
filename = "samples/data1.fq"

def main():
    # read in file and create a list of strings, each string being every
    # 4th line in the file, starting from the second line
    with open(filename, 'r') as file:
        lines = file.readlines()
        lines = lines[1::4]
    # construct a dataframe with columns as the percent A, G, C, T at that specific index
    df = pd.DataFrame(columns=['%A', '%G', '%C', '%T'])
    for i in range(len(lines[0])-1):
        # count the number of each base at that position
        A = sum([1 for line in lines if line[i] == 'A'])
        G = sum([1 for line in lines if line[i] == 'G'])
        C = sum([1 for line in lines if line[i] == 'C'])
        T = sum([1 for line in lines if line[i] == 'T'])
        # calculate the percentage of each base at that position
        total = A + G + C + T
        df.loc[i] = [A/total, G/total, C/total, T/total]
    print(df.head())
    print("Number of rows: ", df.shape[0])
    print(total)

    # plot the percent of each base per read. y axis is percent, and 
    # x axis is read position, or the row numnber
    # write the plot to an html file
    plt.plot(df['%A'], label='A')  
    plt.plot(df['%G'], label='G')
    plt.plot(df['%C'], label='C')
    plt.plot(df['%T'], label='T')
    plt.xlabel('Read Position')
    plt.ylabel('Percent')
    plt.legend()
    plt.show()




if __name__ == "__main__":
    main()