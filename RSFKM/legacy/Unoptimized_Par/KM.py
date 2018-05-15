import numpy as np
import random

MAX_ITERATIONS = 100


def Convergence(OldCentrioids, centroids, iterations):
    if iterations > MAX_ITERATIONS:
        return True
    return np.array_equal(OldCentrioids, centroids)


def GetMembership(DataMatrix, Centroids):
    MinDist = float('Inf')
    MinCentroid = []
    Memberships = []
    SumBuffer = 0.0

    for index, row in enumerate(DataMatrix): #for each object in my data matrix
        for cindex, centroid in enumerate(Centroids): #calculate the distance to each centroid
            #here we get the eucledean distances between our centroid and objects
            for column in range(0, DataMatrix.shape[1]):
                SumBuffer += ( (centroid[column] - row[column]) ** 2 )

            #Square our SumBuffer and store the new minimum distance and Minimum
            # centroid if needed
            #We will use the stored centroid for assignment after all have been compared
            SumBuffer = np.sqrt(SumBuffer)
            if(SumBuffer < MinDist):
                MinDist = SumBuffer
                MinCentroid = cindex

        Memberships.append({"ObjectRow": index, "Centroid": MinCentroid})

        #reset for next row/"object"
        SumBuffer = 0.0
        MinDist = float('Inf')
        MinCentroid = []


    return Memberships


def GetMean(Vectors):
    MeanVector = []
    NumFeatures = Vectors.shape[1]
    NumVectors = Vectors.shape[0]

    for feature in range(0, NumFeatures):
        Sum = 0
        for vector in Vectors:
            Sum += vector[feature]

        MeanVector.append(Sum/NumVectors)

    return MeanVector




def GetCentroids(DataMatrix, Memberships, KClusters):
    Centroids = []
    MeanBuffer = np.empty((0, DataMatrix.shape[1]), float)
    NewCentroid = None

    #load a cluster of associated vectors together into an easily
    # collapasable matrix to find a mean for the centroid
    for centroid in range(0,KClusters):
        for membership in Memberships:
            if(membership["Centroid"] == centroid):
                MeanBuffer = np.vstack([MeanBuffer, DataMatrix[ membership["ObjectRow"] ]])

        #check to ensure that our cetnroid has some objects
        #connected to it
        #Else choose a new random centroid
        if MeanBuffer.shape[0] != 0:
            NewCentroid = GetMean(MeanBuffer)
        else:
            NewCentroid = DataMatrix[random.randint(0,DataMatrix.shape[0])]

        Centroids.append(NewCentroid)

        MeanBuffer = np.empty((0, DataMatrix.shape[1]), float)

    return Centroids


#Initialize KCluster vectors of length numfeatures to act as our centroids
def GetRandomCentroids(DataMatrix, NumFeatures, KClusters):
    Centroids = [];

    for vect in range(0,KClusters):
        Centroids.append( DataMatrix[random.randint(0,DataMatrix.shape[0])] )

    Centroids = np.array(Centroids)

    return Centroids



def KM(DataMatrix, KClusters):
    Centroids = None
    Memberships = []
    DMShape = DataMatrix.shape
    NumFeatures = DMShape[1]

    #bookeeping vars
    Iterations = 0
    OldCentrioids = None

    #Initialize centroids
    Centroids = GetRandomCentroids(DataMatrix, NumFeatures, KClusters)

    #while not convergence
    while not Convergence(OldCentrioids, Centroids, Iterations):
        OldCentrioids = Centroids
        Iterations += 1

        Memberships = GetMembership(DataMatrix, Centroids)
        Centroids = GetCentroids(DataMatrix, Memberships, KClusters)

        print(Iterations)

        #counts number of members in each cluster
        counts = [0] * KClusters

        for Row in Memberships:
            for ndx, Centroid in enumerate(Centroids):
                if( Row["Centroid"] == ndx ):
                    counts[ndx] += 1


        #Should print out the membership counts of each cluster
        for ndx, count in enumerate(counts):
            print("Centroid {}, Members: {}".format(ndx, count))
