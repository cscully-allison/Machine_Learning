import operator
import os
import time
import random
import numpy.linalg as la
import matplotlib.pyplot as plt
import numpy as np
from cvxpy import *



MAX_ITERATIONS = 30 #greater than 30 iterations tends to not give more clear resolution

def RenderMemberships(DataValues, Centroids, MembershipMatrix, Iteration, OutputDirectory):


    x = []
    y = []
    gx = []
    gy = []
    bx = []
    by = []
    wx = []
    wy = []

    #build a color profile
    color = (0.3, 0.0, 0.0)
    bettercolor = (0.3, 0.0, 0.0)
    centroidcolor = (0.0,0.0,0.0)

    for v, Centroid in enumerate(Centroids):
        plt.scatter(Centroid[0], Centroid[1], c=color, marker="x")
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


    if not os.path.isdir(OutputDirectory):
        os.makedirs(OutputDirectory)

    plt.savefig(OutputDirectory + "/{}_Iteration".format(Iteration))
    plt.clf()



def PrintMemberships(Centroids, MembershipMatrix, DataMatrix):
    print "Item, Cluster, Realive Membership [%]"
    for x, item in enumerate(MembershipMatrix):
        for v, membership in enumerate(item):
            if membership > 0.50:
                print  "{0}, {1}, {2}".format(DataMatrix[x], v, membership*100)

    print "Key: "
    for v, centroid in enumerate(Centroids):
        print "Cluster number: ", v, "Centroid: ", centroid, " "

def GetRandomCentroids(DataMatrix, KClusters):
    Centroids = [];
    Selection = [];

    for vect in range(0,KClusters):

        Selection = DataMatrix[random.randint(0,DataMatrix.shape[0]-1)]
        while any((Selection == Centroid).all() for Centroid in Centroids):
            Selection = DataMatrix[random.randint(0,DataMatrix.shape[0]-1)]

        Centroids.append(Selection)

    Centroids = np.array(Centroids)

    return Centroids

def Convergence(OldCentrioids, centroids, iterations):
    if iterations > MAX_ITERATIONS:
        return True
    return np.allclose(OldCentrioids, centroids, rtol=1e-05)


def GetHMatrix(DataMatrix, H, S, V):
    for i, row in enumerate(DataMatrix):
        for k, centroid in enumerate(V):
            H[i][k] = S[i][k] * (la.norm(np.subtract(row, centroid)) ** 2)



def UpdateMembershipMatrix(DataMatrix, H, S, Centroids, MembershipMatrix, RegParam):
    GetHMatrix(DataMatrix, H, S, Centroids)

    #print DataMatrix

    #minimization solving
    #This part essentially corresponds to formula 14
    Ux = Variable(MembershipMatrix.shape[1])
    h_tilde = Parameter(H.shape[1])

    #constraints = [ sum(Ux) == 1 ]
    constraints = [ 0 <= Ux, sum(Ux) == 1]

    for i, Ui in enumerate(MembershipMatrix):
        h_tilde.value = np.multiply(H[i], ( -1/(2*RegParam) ))
        expression = square( norm( Ux - h_tilde ) )
        objective = Minimize(expression)

        prob = Problem(objective, constraints)
        prob.solve()

        if Ux.value is not None:
            #print "Value of Ux", Ux.value
            #print "Value of Ux Transpose", Ux.value.transpose()[0]
            MembershipMatrix[i] = Ux.value.transpose()[0]
            #print "Solution ", prob.value, "Solution Value", Ux.value, "H part", h_tilde.value
        else:
            print "No solution", Ux.value, prob.value, prob.status, h_tilde.value

    #print "This is the U", MembershipMatrix



def UpdateS(DataMatrix, Centroids, S, ThresholdValue):

    for i, row in enumerate(S):
        for k, col in enumerate(row):
            if la.norm(np.subtract(DataMatrix[i], Centroids[k])) == 0:
                S[i][k]
            elif la.norm(np.subtract(DataMatrix[i], Centroids[k])) > ThresholdValue:
                S[i][k] = 0
            else:
                S[i][k] = 1/( la.norm(np.subtract(DataMatrix[i], Centroids[k])) )

    #print "This is the S", S




# v = centroids, s= is matrix holding s aux vars, and u is membership matrix
def FindCentroids(DataMatrix, V, S, U):
    SummedDenom = 0.0 #this is the buffer for the denomonator of our vector centroid function
    Scalar = 0.0

    #for each centroid
    for k, vk in enumerate(V):
        tempCentroid = np.zeros([DataMatrix.shape[1]])

        #for each item in the data matrix
        for i, row in enumerate(DataMatrix):
            Scalar = S[i][k] * U[i][k] #get the scalar we are multiplying our x_i by

            if Scalar != 0.0:

                #multiply our scalar against data matrix vector i
                # and add the result to temp centroid f
                #tempCentroid = np.multiply(DataMatrix[i], Scalar)
                for f, feature in enumerate(vk):
                    tempCentroid[f] += Scalar * DataMatrix[i][f]

            #keep track of our denomonator sum
            SummedDenom += Scalar


        #divide our centroid vector by the denomonator we computed
        if SummedDenom != 0.0:
            tempCentroid = np.multiply(tempCentroid, 1/SummedDenom)


        #iterate over number of features in our Vectors jsut to copy one to the other
        V[k] = np.copy(tempCentroid)

        #reset our temp centroid
        SummedDenom = 0.0

    return V




# Function: RSFKM
# Description: Robust and Sparse Fuzzy K-Means Algorithim (to be implemented and modified for data imputation)
# Parameter 1: DataMatrix [NumPy Array] : An array contianing all my Data values (Rows are items, cols are features)
# Parameter 2: KClusters [int] : The number of clusters we are grouping into
# Parameter 3: RegParam [float] : A regularization parameter that puts a restriction on the minimum distance between a data point
#                                  and a cluster's center and prevents membership from having extreme values, 0 and 1
# Parameter 4: ThresholdValue [float] : Controls the number of outliers.  If the residual of a sample to centroid is larger than <ThresholdValue>,
#                                       it is re-garded as outlier and not used to learn centroid matrix V since the corresponding s_ik is zero"
def RSFKM(DataMatrix, KClusters, RegParam, ThresholdValue, OutputDirectory):
    #variables
    Centroids = np.empty([KClusters, DataMatrix.shape[1]], dtype=float) #corresponds to V in paper
    MembershipMatrix = np.empty([DataMatrix.shape[0], KClusters], dtype=float) #should be of shape i rows and k columns; corresponds to U in paper
    MatrixH = np.empty([DataMatrix.shape[0], KClusters], dtype=float) #derived from our s aux variable
    S = np.empty([DataMatrix.shape[0], KClusters], dtype=float) #holds our calulated s_ik
    TimeStep = 0
    OldCentrioids = np.empty([KClusters, DataMatrix.shape[1]], dtype=float)

    #initalization
    for rndx, row in enumerate(DataMatrix):
        for col in range(0, KClusters):
            MemebershipMatrix = 1/KClusters
            # if col == 0:
            #     MembershipMatrix[rndx][col] = 0.9
            # elif col == 1:
            #     MembershipMatrix[rndx][col] = 0.08
            # else:
            #     MembershipMatrix[rndx][col] = 0.0001 #this needs to be fixed to enforce constraint of U1 = 1 in cases of odd number of clusters like 3

            S[rndx][col] = 1

    #Centroids = FindCentroids(DataMatrix, Centroids, S, MembershipMatrix)
    Centroids = GetRandomCentroids(DataMatrix, KClusters)
    print Centroids



    while not Convergence(OldCentrioids, Centroids, TimeStep):

        OldCentrioids = np.copy(Centroids)

        UpdateMembershipMatrix(DataMatrix, MatrixH, S, Centroids, MembershipMatrix, RegParam)
        Centroids = FindCentroids(DataMatrix, Centroids, S, MembershipMatrix)
        UpdateS(DataMatrix, Centroids, S, ThresholdValue)

        print "Iteration ", TimeStep, " complete"
        TimeStep += 1

        RenderMemberships(DataMatrix, Centroids, MembershipMatrix, TimeStep, OutputDirectory)

    return { "U": MembershipMatrix, "V":Centroids  }

def ImputeData(DataMatrix, KClusters, RegParam, ThresholdValue):
    comparision_vect = []

    #strip away the first feature from n rows and store for later
    for row in DataMatrix:
        comparision_vect.append(row[0])

    #perform clustering using all data but only using n-1 features

    #take returned U and V and find the geo average for feature 1 in cluster v^i and use that as our imputed Data

    #compare accuracy of imputed value and real value
    return
