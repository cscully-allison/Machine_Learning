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
    vx = []
    vy = []

    Data = ReadInData("smol2fd.csv")
    DataValues = Data["DataFrame"]

    CleanData(DataValues)

    #KM(DataValues,20)

    UVBundle = RSFKM(DataValues, 4, 3, 10)
    MembershipMatrix = UVBundle["U"]
    Centroids = UVBundle["V"]

    #print MembershipMatrix
    print Centroids

    #build a color profile
    color = (0.5, 0.0, 0.0)
    bettercolor = (0.6, 0.0, 0.0)
    centroidcolor = (0.0,0.0,0.0)

    for v, Centroid in enumerate(Centroids):
        plt.scatter(Centroid[0], Centroid[1], c=centroidcolor, marker="x")
        centroidcolor = (centroidcolor[0]+.10, centroidcolor[1]+.10,centroidcolor[2]+.10)
        for xi, item in enumerate(DataValues):
            if MembershipMatrix[xi][v] > 0.9:
                vx.append(item[0])
                vy.append(item[1])
            elif MembershipMatrix[xi][v] > 0.5:
                x.append(item[0])
                y.append(item[1])
        plt.scatter(vx,vy, c=bettercolor, marker="o")
        plt.scatter(x, y, c=color, marker=".")
        x = vx = []
        y = vy = []


        color = (color[2]+.1,color[0],color[1]+.1)
        bettercolor = (bettercolor[2]+.1,bettercolor[0],bettercolor[1]+.1)



    plt.show()




main();
