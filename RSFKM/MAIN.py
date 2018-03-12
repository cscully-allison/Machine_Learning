import csv
import numpy as np
import random
from KM import KM
from RSFKM import RSFKM
import matplotlib.pyplot as plt


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

def CleanData(DataValues):
    for row in DataValues:
        for col in row:
            if col is None:
                print "ajksff"



def main():
    Data = None
    DataValues = None
    UVBundle = None
    MembershipMatrix = None
    Centroids = None
    x = []
    y = []

    Data = ReadInData("smol2fd.csv")
    DataValues = Data["DataFrame"]

    CleanData(DataValues)

    #KM(DataValues,20)

    UVBundle = RSFKM(DataValues, 2, 3, 10)
    MembershipMatrix = UVBundle["U"]
    Centroids = UVBundle["V"]

    for item in DataValues:
        x.append(item[0])
        y.append(item[1])

    for centroid in DataValues:
        x.append(centroid[0])
        y.append(centroid[1])


    plt.plot(x, y)



main();
