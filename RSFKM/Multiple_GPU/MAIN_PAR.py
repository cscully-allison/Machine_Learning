#!/usr/bin/env python
import time
import operator
import sys
import Queue

from threading import Thread

import numpy as np
import pycuda.driver as drv
from RSFKM_PAR import threadedRSFKM, ReadInData


"""Collect command-line options in a dictionary"""

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts




def main():

    print "get here?\n"

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
    MissingPercent = float(Args["-P"])

    queue = Queue.Queue()
    Threads = []
    Returns = []

    #
    # #We need to call this on n threads and pass in a device ID
    print drv.Device.count()
    #
    for DID in range(drv.Device.count()):
    #     #call our threaded function and pass all relevant data as a "struct"
          Threads.append(Thread( target=threadedRSFKM, args=({"DID": DID, "DataSource":DataSource, "DTColFlag":DTColFlag, "NumRows":NumRows, "NumCols":NumCols, "NumClusters": NumClusters, "RegParam": RegParam, "ThresholdValue":ThresholdValue, "OutputDirectory": OutputDirectory, "MissingPercent": MissingPercent} , queue) ))
    #
    for thread in Threads:
        thread.start()
    #
    for thread in Threads:
        thread.join()
    #
    while not queue.empty():
         Returns.append(queue.get())
    #
    for gpu, data in enumerate(Returns):
        print "GPU:", data["DID"], data["RSME"], data["Output"]
        # for i, value in enumerate(data["GuessedValues"]):
        #     print "Guess:", data["GuessedValues"][i], " vs. Real: ", data["GroundTruths"][i]
    #
    # print "Happens last right!?"







    #Main Driver of RSFKM (ROBUST AND SPARSE FUZZY K MEANS)
    # UVBundle = ImputeData(DataValues, NumClusters, RegParam, ThresholdValue, OutputDirectory, MissingPercent)



    # MembershipMatrix = UVBundle["U"]
    # Centroids = UVBundle["V"]

    #print MembershipMatrix
    #print Centroids
    #PrintMemberships(Centroids,MembershipMatrix,DataValues)
    #RenderMemberships(DataValues, Centroids, MembershipMatrix, 0, OutputDirectory)




main();
