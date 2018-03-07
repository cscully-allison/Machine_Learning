import csv
import numpy as np
import random
from KM import KM
from RSFKM import RSFKM


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

    Data = ReadInData("smoldata.csv")
    DataValues = Data["DataFrame"]


    #KM(DataValues,20)
    RSFKM(DataValues, 5, 0.5, 3)


main();
