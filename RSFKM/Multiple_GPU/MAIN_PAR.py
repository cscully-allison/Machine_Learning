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
    Data = None
    DataValues = None
    data = None
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
    TimeAVG = 0.0
    TimePerIterAVG = 0.0
    IterAVG = 0
    RSMEAVG = 0.0

    #
    # #We need to call this on n threads and pass in a device ID
    for DID in range(drv.Device.count()):
        #call our threaded function and pass all relevant data as a "struct"
        Threads.append(Thread( target=threadedRSFKM, args=({"DID": DID, "DataSource":DataSource, "DTColFlag":DTColFlag, "NumRows":NumRows, "NumCols":NumCols, "NumClusters": NumClusters, "RegParam": RegParam, "ThresholdValue":ThresholdValue, "OutputDirectory": OutputDirectory, "MissingPercent": MissingPercent} , queue) ))

    for thread in Threads:
        thread.start()

    for thread in Threads:
        thread.join()

    while not queue.empty():
         Returns.append(queue.get())

    for gpu, data in enumerate(Returns):
        TimeAVG += data["Time(ms)"]
        TimePerIterAVG += data["Time(ms)"]/data["Iter"]
        IterAVG += data["Iter"]
        RSMEAVG += data["RSME"]

        # print "GPU:", data["DID"], data["RSME"]
        # for i, value in enumerate(data["GuessedValues"]):
        #     print "Guess:", data["GuessedValues"][i], " vs. Real: ", data["GroundTruths"][i]

    TimeAVG = TimeAVG/drv.Device.count()
    TimePerIterAVG = TimePerIterAVG/drv.Device.count()
    IterAVG = IterAVG/drv.Device.count()
    RSMEAVG = RSMEAVG/drv.Device.count()

    print "{},{},{},{},{},{},{},{},{}".format(drv.Device.count(), IterAVG, NumRows, NumRows*drv.Device.count(), NumCols, NumClusters,RSMEAVG, TimeAVG, TimePerIterAVG)






main();
