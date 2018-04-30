import operator
import os
import time
import random
import math
import numpy.linalg as la
#import matplotlib.pyplot as plt
import numpy as np
import pycuda.autoinit as ai
import pycuda.driver as drv
from pycuda.compiler import SourceModule
from cvxpy import *



MAX_ITERATIONS = 30 #greater than 30 iterations tends to not give more clear resolution

def RenderMemberships(DataValues, Centroids, MembershipMatrix, Iteration, OutputDirectory):
    # y = []
    # x = []
    # gx = []
    # gy = []
    # bx = []
    # by = []
    # wx = []
    # wy = []
    #
    # #build a color profile
    # color = (0.3, 0.0, 0.0)
    # bettercolor = (0.3, 0.0, 0.0)
    # centroidcolor = (0.0,0.0,0.0)
    #
    # for v, Centroid in enumerate(Centroids):
    #     plt.scatter(Centroid[0], Centroid[1], c=color, marker="x")
    #     for ui, item in enumerate(MembershipMatrix):
    #         index, membership = max(enumerate(item), key=operator.itemgetter(1))
    #         #print membership
    #         if index == v:
    #             if membership > 0.9:
    #                 gx.append(DataValues[ui][0])
    #                 gy.append(DataValues[ui][1])
    #             elif membership > 0.5:
    #                 x.append(DataValues[ui][0])
    #                 y.append(DataValues[ui][1])
    #             elif membership > 0.2:
    #                 bx.append(DataValues[ui][0])
    #                 by.append(DataValues[ui][1])
    #             elif membership >= 0.0:
    #                 wx.append(DataValues[ui][0])
    #                 wy.append(DataValues[ui][1])
    #     plt.scatter(gx,gy, c=bettercolor, marker="o")
    #     plt.scatter(x, y, c=color, marker=".")
    #     plt.scatter(bx, by, c=color, marker=",")
    #     plt.scatter(wx, wy, c=(0.0,0.0,0.0), marker="h")
    #
    #     x = []
    #     y = []
    #     gx = []
    #     gy = []
    #     bx = []
    #     by = []
    #     wx = []
    #     wy = []
    #
    #
    #     color = (color[2]+.05,color[0],color[1]+.05)
    #     bettercolor = (bettercolor[2]+.05,bettercolor[0],bettercolor[1]+.05)
    #
    #
    # if not os.path.isdir(OutputDirectory):
    #     os.makedirs(OutputDirectory)
    #
    # plt.savefig(OutputDirectory + "/{}_Iteration".format(Iteration))
    # plt.clf()
    return



def PrintMemberships(Centroids, MembershipMatrix, DataMatrix):
    print "Item, Cluster, Realive Membership [%]"
    for x, item in enumerate(MembershipMatrix):
        for v, membership in enumerate(item):
            if membership > 0.50:
                print  "{0}, {1}, {2}".format(DataMatrix[x], v, membership*100)

    print "Key: "
    for v, centroid in enumerate(Centroids):
        print "Cluster number: ", v, "Centroid: ", centroid, " "
        count = 0
        for item in MembershipMatrix:
            if item[v] > .50:
                count += 1
        print "Number of items in cluster: ", count


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


#@profile
def GetHMatrix(DataMatrix, H, S, V):
    for i, row in enumerate(DataMatrix):
        for k, centroid in enumerate(V):
            H[i][k] = S[i][k] * ( la.norm(np.subtract(row, centroid)) ** 2)


#@profile
def UpdateMembershipMatrix(DataMatrix_GPU, DataMatrix, H_GPU, H, S_GPU, S, Centroids_GPU, V, MembershipMatrix_GPU, MembershipMatrix, RegParam, mod):

    build_h_matrix = mod.get_function("build_h_matrix");
    call_solver = mod.get_function("update_membership_matrix");


    #Declare some temp variable
    H_Flat = H.flatten().astype(np.float64)
    MM_Flat = MembershipMatrix.flatten().astype(np.float64)


    NumRows = np.int32(H.shape[0])
    NumCentroids = np.int32(H.shape[1])
    Int32NumFeatures = np.int32(DataMatrix.shape[1])


    StrictRegParam = np.float64(RegParam)



    #GetHMatrix(DataMatrix, H, S, V)



    """Builds up H Matrix on the Cuda Device"""
    build_h_matrix(H_GPU, DataMatrix_GPU, S_GPU, Centroids_GPU, block=(DataMatrix.shape[1],1,1), grid=(H.shape[1], H.shape[0], 1))


    ai.context.synchronize()

    drv.memcpy_dtoh(H_Flat, H_GPU)

    H_Flat = np.reshape(H_Flat, H.shape)



    """Do some quick maths to determine the numer of blocks needed"""
    """We will oprimize this part in the future because it causes problems
    with  warp effiency, There are a dozen+ warps which are running uselessly :/"""

    call_solver(MembershipMatrix_GPU, H_GPU, StrictRegParam, NumCentroids, NumRows, block=(1,1,1), grid=(DataMatrix.shape[0],1))
    # call_solver(drv.Out(MM_Flat), MembershipMatrix_GPU, H_GPU, StrictRegParam, NumCentroids, NumRows, block=(1,1,1), grid=(DataMatrix.shape[0],1)) #TPB cannot exceed 256

    ai.context.synchronize()


    # MM_Flat = np.reshape(MM_Flat, MembershipMatrix.shape)
    #
    #
    # #doing this to ensure that the membership matrix gets copied out
    # for i,row in enumerate(MembershipMatrix):
    #      MembershipMatrix[i] = MM_Flat[i]
    #
    # for i in MembershipMatrix:
    #     print sum(i)
    #


#@profile
def UpdateS(DataMatrix, Centroids, S, ThresholdValue):
    NormResult = None

    for i, row in enumerate(S):
        for k, col in enumerate(row):
            NormResult = la.norm(np.subtract(DataMatrix[i], Centroids[k]))
            if NormResult == 0:
                S[i][k]
            elif NormResult > ThresholdValue:
                S[i][k] = 0
            else:
                S[i][k] = 1/NormResult





# v = centroids, s= is matrix holding s aux vars, and u is membership matrix
#@profile
def FindCentroids( DataMatrix, V, S, U):
    SummedDenom = 0.0 #this is the buffer for the denomonator of our vector centroid function
    Scalar = 0.0

    #for each centroid
    for k, vk in enumerate(V):
        tempCentroid = np.zeros([DataMatrix.shape[1]])

        #for each item in the data matrix
        for i, row in enumerate(DataMatrix):
            Scalar = S[i][k] * U[i][k] #get the scalar we are multiplying our x_i by


            #if Scalar != 0.0:

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

        # if k == 3:
        #     print tempCentroid, "Summed Denom: ", SummedDenom, 1/SummedDenom


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
    MatrixH = np.zeros([DataMatrix.shape[0], KClusters], dtype=float) #derived from our s aux variable
    S = np.zeros([DataMatrix.shape[0], KClusters], dtype=float) #holds our calulated s_ik
    OldCentrioids = np.empty([KClusters, DataMatrix.shape[1]], dtype=float)


    TimeStep = 0

    #allocate some space on the kernel
    DataMatrix_GPU = drv.mem_alloc(DataMatrix.flatten().astype(np.float64).nbytes)
    MembershipMatrix_GPU = drv.mem_alloc(MembershipMatrix.flatten().astype(np.float64).nbytes)
    Centroids_GPU = drv.mem_alloc(Centroids.flatten().astype(np.float64).nbytes)
    S_GPU = drv.mem_alloc(S.flatten().astype(np.float64).nbytes)
    AUX_GPU = drv.mem_alloc(MatrixH.flatten().astype(np.float64).nbytes)

    THREADS = 640



    drv.memcpy_htod(DataMatrix_GPU, DataMatrix.flatten().astype(np.float64))

    #create module 1 time
    source = open("cuda_solver/solver.cu",'r')
    sourceCode = source.read()

    mod = SourceModule(sourceCode)


    init_S = mod.get_function("init_S");
    load_scalar_buffer = mod.get_function("load_scalar_buffer")
    calculate_centroids = mod.get_function("calculate_centroids")
    find_centroids = mod.get_function("find_centroids")
    update_S = mod.get_function("update_S")





    #initalization
    for rndx, row in enumerate(DataMatrix):
        for col in range(0, KClusters):
            #MemebershipMatrix = 1/KClusters
            if col == 0:
                MembershipMatrix[rndx][col] = 0.9
            elif col == 1:
                MembershipMatrix[rndx][col] = 0.08
            else:
                MembershipMatrix[rndx][col] = 0.0001 #this needs to be fixed to enforce constraint of U1 = 1 in cases of odd number of clusters like 3

            S[rndx][col] = 1.0




    # initialize s on the gpu
    # may not even be worth it
    # S.shape[0] is number of rows
    #init_S(S_GPU, np.int32(KClusters), block=(KClusters,1,1), grid=(S.shape[0],1,1));

    drv.memcpy_htod(S_GPU, S.flatten().astype(np.float64))

    ai.context.synchronize()

    Centroids = GetRandomCentroids(DataMatrix, KClusters)

    drv.memcpy_htod(MembershipMatrix_GPU, MembershipMatrix.flatten().astype(np.float64))
    drv.memcpy_htod(Centroids_GPU, Centroids.flatten().astype(np.float64))



#CORE PROCESSING LOOP! ----------------------------------

    while not Convergence(OldCentrioids, Centroids, TimeStep):
        OldCentrioids = np.copy(Centroids)

        UpdateMembershipMatrix(DataMatrix_GPU, DataMatrix, AUX_GPU, MatrixH, S_GPU, S, Centroids_GPU, Centroids,  MembershipMatrix_GPU, MembershipMatrix, RegParam, mod)


        #going to use H matrix as a scalar buffer to save on space
        load_scalar_buffer(AUX_GPU, S_GPU, MembershipMatrix_GPU, np.int32(DataMatrix.shape[0]), np.int32(Centroids.shape[0]), block=(THREADS, 1, 1), grid=(Centroids.shape[0],1,1))
        ai.context.synchronize()
        calculate_centroids(DataMatrix_GPU, Centroids_GPU, AUX_GPU, np.int32(DataMatrix.shape[0]), np.int32(DataMatrix.shape[1]), block=(THREADS, 1, 1), grid=(Centroids.shape[0], Centroids.shape[1]), shared=(DataMatrix.shape[0] * 8))
        ai.context.synchronize()



        find_centroids(DataMatrix_GPU, Centroids_GPU, AUX_GPU, np.int32(DataMatrix.shape[0]), np.int32(DataMatrix.shape[1]), block=(32, 1, 1), grid=( Centroids.shape[0], 1))
        ai.context.synchronize()


        #Centroids = FindCentroids(DataMatrix, Centroids, S, MembershipMatrix)

        update_S( S_GPU, DataMatrix_GPU, Centroids_GPU, np.float64(ThresholdValue), np.int32(DataMatrix.shape[0]), np.int32(DataMatrix.shape[1]), block=(DataMatrix.shape[1],1,1), grid=(DataMatrix.shape[0], Centroids.shape[0], 1), shared=(DataMatrix.shape[1] * 8))
        ai.context.synchronize()

        #UpdateS(DataMatrix, Centroids, S, ThresholdValue)

        Centroids_Host = Centroids.flatten().astype(np.float)
        drv.memcpy_dtoh(Centroids_Host, Centroids_GPU)
        Centroids = np.reshape(Centroids_Host, Centroids.shape)


        print TimeStep

        TimeStep += 1



    MembershipMatrix = MembershipMatrix.flatten().astype(np.float64)
    drv.memcpy_dtoh(MembershipMatrix, MembershipMatrix_GPU)
    MembershipMatrix = np.reshape(MembershipMatrix, (DataMatrix.shape[0], Centroids.shape[0]))




    #PrintMemberships(Centroids, MembershipMatrix, DataMatrix)
#        RenderMemberships(DataMatrix, Centroids, MembershipMatrix, TimeStep, OutputDirectory)

    return { "U": MembershipMatrix, "V":Centroids, "Iter":TimeStep  }





def ImputeData(DataMatrix, KClusters, RegParam, ThresholdValue, OutputDirectory, MissingPercent):
    NumRows = DataMatrix.shape[0]
    NumCols = DataMatrix.shape[1]

    NumMissing = int(DataMatrix.shape[0]*MissingPercent)

    ComparisionVect = []
    ExistingValues = []
    ClusterMapping = np.zeros((DataMatrix.shape[0]), dtype=int)
    ImputableDataValues = []

    for i in range(KClusters):
        ImputableDataValues.append({"NumValues":0,"Sum":0.0})

    #strip away the first feature from n rows and store for later
    for row in DataMatrix:
        ComparisionVect.append(row[0])

    for i in range(0, int(DataMatrix.shape[0])):
        if(i > NumMissing):
            ExistingValues.append(ComparisionVect[i])
        else:
            ExistingValues.append(float('nan'))


    #perform clustering using all data but only using n-1 features
    DataMatrix = np.delete(DataMatrix, 0, axis=1)

    start = time.time()
    UVBundle = RSFKM(DataMatrix, KClusters, RegParam, ThresholdValue, OutputDirectory)
    end = time.time()



    #take returned U and V and find the geo average for feature 1 in cluster v^i and use that as our imputed Data
    #Load 15 clusters with sums of exising data
    for i, line in enumerate(UVBundle["U"]):
        for k, centroid in enumerate(UVBundle["V"]):
            if line[k] > .5:
                ClusterMapping[i] = k
            if line[k] > .5 and not math.isnan(ExistingValues[i]):
                ImputableDataValues[k]["NumValues"] += 1
                ImputableDataValues[k]["Sum"] += ExistingValues[i]

    for item in ImputableDataValues:
        if item["NumValues"] > 0:
            item["Avg"] = item["Sum"]/item["NumValues"]

    for value in ImputableDataValues:
        print value

    #impute data
    for i, value in enumerate(ExistingValues):
        if math.isnan(value):
            #find imputable data point from imputable data values (derefrenced via cluster mapping)
            ExistingValues[i] = ImputableDataValues[ClusterMapping[i]]["Avg"]
            #print "Guess:", ExistingValues[i], " vs. Real: ", ComparisionVect[i]


    sum = 0
    for i in range(NumMissing):
        sum += ((ExistingValues[i] - ComparisionVect[i]) ** 2)

    RSME = math.sqrt(sum/NumMissing)

    print "RSME: ", RSME
    print "{},{},{},{},{},{}".format( UVBundle["Iter"], NumRows, NumCols, KClusters, ((end - start)*1000), ((end - start)*1000)/UVBundle["Iter"])

    return UVBundle;
