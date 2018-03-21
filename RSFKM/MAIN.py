import operator
import csv
import numpy as np
import random
import sys
from KM import KM
from RSFKM import RSFKM, RenderMemberships

"""Collect command-line options in a dictionary"""

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts


def ReadInData(source):
    TextValues = [] #row major matrix representing the set of all data points we will be grouping
    FloatValues = []
    DtMappings = [] #mappings of date time values to thier associated row of data in the object set

    with open(source, 'rt') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if '' not in row: #throw out rows with missing data
                DtMappings.append(row[0])
                TextValues.append(row[1:len(row)])


    TextValues = np.array(TextValues)
    FloatValues = TextValues.astype(np.float)

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
    RegParam = int(Args["-r"])
    ThresholdValue = int(Args["-t"])

    print Args

    Data = ReadInData(DataSource)
    DataValues = Data["DataFrame"]

    #CleanData(DataValues)

    #KM(DataValues,20)

    #UVBundle = RSFKM(DataValues, 15, 8, 20, OutputDirectory)
    UVBundle = RSFKM(DataValues, NumClusters, RegParam, ThresholdValue, OutputDirectory)
    MembershipMatrix = UVBundle["U"]
    Centroids = UVBundle["V"]

    #print MembershipMatrix
    print Centroids
    RenderMemberships(DataValues, Centroids, MembershipMatrix, 0, OutputDirectory)




main();
