#!/usr/bin/env python
import time
import operator
import csv
import numpy as np
import random
import sys
from KM import KM
from RSFKM_PAR import RSFKM, RenderMemberships, PrintMemberships


"""Collect command-line options in a dictionary"""

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts


def ReadInData(source, DTColFlag, NumRows, NumCols):
    TextValues = [] #row major matrix representing the set of all data points we will be grouping
    FloatValues = []
    DtMappings = [] #mappings of date time values to thier associated row of data in the object set

    with open(source, 'rt') as csvfile:
        reader = csv.reader(csvfile)
        for i, row in enumerate(reader):
            if i >= NumRows:
                break
            if len(row) < NumCols+1:
                NumCols = len(row)
            if '' not in row: #throw out rows with missing data
                if DTColFlag is 1:
                    DtMappings.append(row[0])
                    TextValues.append(row[1:NumCols+1])
                else:
                    TextValues.append(row)



    TextValues = np.array(TextValues)
    FloatValues = TextValues.astype(np.float32)

    return {"DTMappings": DtMappings, "DataFrame": FloatValues}




def main():
    Data = None
    DataValues = None
    UVBundle = None
    MembershipMatrix = None
    Centroids = None

    #retrieve passed in args
    Args = getopts(sys.argv)

    DataSource = Args["-i"]
    OutputDirectory = Args["-o"]
    NumClusters = int(Args["-k"])
    RegParam = float(Args["-r"])
    ThresholdValue = float(Args["-t"])
    DTColFlag = int(Args["-l"])
    NumRows = int(Args["-rw"])
    NumCols =  int(Args["-c"])

    Data = ReadInData(DataSource, DTColFlag, NumRows, NumCols)
    DataValues = Data["DataFrame"]

    #UVBundle = RSFKM(DataValues, 15, 8, 20, OutputDirectory)

    #Main Driver of RSFKM (ROBUST AND SPARSE FUZZY K MEANS)

    start = time.time()
    UVBundle = RSFKM(DataValues, NumClusters, RegParam, ThresholdValue, OutputDirectory)
    end = time.time()

    print "{},{},{},{},{},{}".format( UVBundle["Iter"], NumRows, NumCols, NumClusters, ((end - start)*1000), ((end - start)*1000)/UVBundle["Iter"] )


    MembershipMatrix = UVBundle["U"]
    Centroids = UVBundle["V"]

    #print MembershipMatrix
    #print Centroids
    #PrintMemberships(Centroids,MembershipMatrix,DataValues)
    #RenderMemberships(DataValues, Centroids, MembershipMatrix, 0, OutputDirectory)




main();
