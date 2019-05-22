#!/usr/bin/python

import sys
import numpy as np
import os
import operator
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys, getopt

## Argument parser

import argparse

parser = argparse.ArgumentParser(description='Using the count matrix generates venn diagram with genes that have mean normalixzed counts > 10 for each group.')
parser.add_argument('-i', metavar='<inputfile>', help='input count file', required=True)
parser.add_argument('-m', metavar='<metadatafile>', help='metafile with description of samples', required=True)
parser.add_argument('-s', metavar='<sample>', help='name given to sample file', required=True)
parser.add_argument('-p', metavar='<pattern>', help='name given to pattern file', required=True)
parser.add_argument('-o', metavar='<outfiel>', help='figure name', required=True)

args = parser.parse_args()


inputfile = args.i 
metadatafile = args.m 
sample = args.s 
pattern = args.p 
outfile = args.o


# Metadta file with all samples information

def metadata(meta):
    header = []
    with open(meta) as f:
        for line in f:
            samp = line.strip()
            header.append(samp)
    return header

metaD = metadata(metadatafile)

# get the index of samples

def index():
    d = {}
    with open(inputfile) as f:
        header = f.readline().strip().split('\t')
        for i in metaD:
            d[i] = [header.index(s) for s in header if i in s]
    return d

# Get the genes that have mean normalized value > 10 for all the smaples.

def counts(sample):
    avg = []
    genelst = []
    countAvg = 0
    with open(inputfile) as f:
        header = f.readline()
        lines = f.readlines()
        
        for line in lines:
            columns = line.split('\t')
            for vv in index()[sample]:
                val = np.log2(int(columns[vv])+0.0001)
                avg.append(float(val))
            mean = np.mean(avg)
            if mean > 10.0:
                countAvg +=1
                genelst.append(columns[0])
            #print avg
            avg =[]
    return genelst

# Generate a file with genes in different samples
def write():
    import csv
    import itertools as it
    biglst = []
    with open(sample, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(metaD)
        for i, val in enumerate(metaD):
            biglst.append(counts(val))
            csv.writer(f).writerows(it.izip_longest(*biglst))

write()

def geneNames():
    geneLst = []
    with open(sample) as f:
        header = f.readline().strip(',')
        for line in f:
            ln = line.strip().split(',')
            #print ln
            geneLst.append(ln)

    geneLst1 = [item for sublist in geneLst for item in sublist]
    newLst = list(set(geneLst1))
    return newLst


def getIndex(gn):
    with open(sample) as f1:
        allLst = []
        header = f1.readline()
        lines = f1.readlines()
        for line1 in lines:
            if gn in line1: #and gn =='':
                ln1 = line1.strip().split(',')
                #print ln1
                #mock = [0 for j in range(5)]
                new = [ind for ind, val in enumerate(ln1) if val==gn]
                allLst.append(new)
    
    return allLst

# generate pattern file
def countPattern():
    d1 = {}
    final = []
    finalD = {}
    myfile=open(pattern,'w')
    myfile.write('\t'+'\t'.join(metaD) + '\t' + 'Weight\n')
    for gene in geneNames():
        ll = [item for sublist in getIndex(gene) for item in sublist]
        mock = [0 for i in range(5)]
        for j in ll:
            mock[j] = 1
        mm = [str(jj) for jj in mock]
        strF = ''.join(mm)
        final.append(strF)

    for m in final:
       finalD[m] = finalD.get(m, 0) + 1

    for k, v in sorted(finalD.iteritems()):
        finalS = '\t'.join(list(k))
        myfile.write(str(k) + '\t' + finalS +'\t' + str(v) + '\n')
    myfile.close()

    
countPattern()





######################################
#for scatter plot
######################################

def plot():
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np

    import matplotlib.gridspec as gridspec


    with open(pattern) as f:
        header = f.readline().strip().split('\t')
        lines = f.readlines()
        lst = []
        wt = []
        newlst = []
        for line in lines:
            ln = line.strip().split('\t')
            wt.append(int(ln[-1]))
            ln1 = [int(i) for i in ln[1:]]
            lst.append(ln1[:-1])
   
            for j in lst:
                newlst.append([i for i, s in enumerate(j) if s==1])
        # set the spacing between axes
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.figure(1)
        plt.subplot(211)
        centersScatter = np.arange(len(wt))
        width = 1
        plt.bar(centersScatter + width, wt, width, color = 'red', edgecolor = "k")
        plt.xticks([])


        plt.figure(1)
        plt.subplot(212)

    ######################################
    #for Bar plot
    ######################################

    centersBar = np.arange(len(newlst))

    for i, j in zip(centersBar, newlst[:-1]):
        for k in range(len(j)):
            plt.plot(i, j[k], 'o', color='red', markersize=4, markeredgecolor='k', linewidth=0.01)

    centers = np.arange(len(header[:-1]))
    plt.yticks(centers, header[:-1])
    
    plt.xticks(centersBar, rotation = 90, fontsize=7)
    

    plt.xticks([])
    plt.grid()
    plt.show()
    plt.savefig(outfile+'.png')

plot()


