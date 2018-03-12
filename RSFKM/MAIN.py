import operator
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
    # if MembershipMatrix[xi][v] > 0.9:
            #     vx.append(item[0])
            #     vy.append(item[1])
            # elif MembershipMatrix[xi][v] > 0.5:
            #     x.append(item[0])
            #     y.append(item[1])
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
    gx = []
    gy = []
    bx = []
    by = []
    wx = []
    wy = []

    Data = ReadInData("smol2fd.csv")
    DataValues = Data["DataFrame"]

    CleanData(DataValues)

    #KM(DataValues,20)

    UVBundle = RSFKM(DataValues, 10, 3, 10)
    MembershipMatrix = UVBundle["U"]
    Centroids = UVBundle["V"]

    #print MembershipMatrix
    print Centroids

    #build a color profile
    color = (0.3, 0.0, 0.0)
    bettercolor = (0.3, 0.0, 0.0)
    centroidcolor = (0.0,0.0,0.0)

    for v, Centroid in enumerate(Centroids):
        plt.scatter(Centroid[0], Centroid[1], c=centroidcolor, marker="x")
        for ui, item in enumerate(MembershipMatrix):
            index, membership = max(enumerate(item), key=operator.itemgetter(1))
            #print membership
            if index == v:
                if membership > 0.9:
                    gx.append(DataValues[ui][0])
                    gy.append(DataValues[ui][1])
                elif membership > 0.5:
                    x.append(DataValues[ui][0])
                    y.append(DataValues[ui][1])
                elif membership > 0.2:
                    bx.append(DataValues[ui][0])
                    by.append(DataValues[ui][1])
                elif membership >= 0.0:
                    wx.append(DataValues[ui][0])
                    wy.append(DataValues[ui][1])
        plt.scatter(gx,gy, c=bettercolor, marker="o")
        plt.scatter(x, y, c=color, marker=".")
        plt.scatter(bx, by, c=color, marker=",")
        plt.scatter(wx, wy, c=(0.0,0.0,0.0), marker="h")

        x = []
        y = []
        gx = []
        gy = []
        bx = []
        by = []
        wx = []
        wy = []


        color = (color[2]+.05,color[0],color[1]+.05)
        bettercolor = (bettercolor[2]+.05,bettercolor[0],bettercolor[1]+.05)



    plt.show()




main();
