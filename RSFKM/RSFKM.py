import numpy as np




# v = centroids, s= is matrix holding s aux vars, and u is membership matrix
def FindCentroids(DataMatrix, V, S, U):
    tempCentroid = np.zeros([DataMatrix.shape[1]])
    SummedDenom = 0.0 #this is the buffer for the denomonator of our vector centroid function
    Scalar = 0.0


    #for each centroid
    for k, vk in enumerate(V):
        #for each item in the data matrix
        for i, row in enumerate(DataMatrix):
            Scalar = S[i][k] * U[i][k] #get the scalar we are multiplying our x_i by

            #multiply our scalar against data matrix vector i
            # and add the result to temp centroid f
            for f, feature in enumerate(vk):
                tempCentroid[f] += Scalar * DataMatrix[i][f]

            #keep track of our denomonator sum
            SummedDenom += Scalar

        #divide our centroid vector by the denomonator we computed
        tempCentroid = np.multiply(tempCentroid, 1/SummedDenom)


        #iterate over number of features in our Vectors jsut to copy one to the other
        V[k] = np.copy(tempCentroid)

        #reset our temp centroid
        tempCentroid = np.zeros([DataMatrix.shape[1]])
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
def RSFKM(DataMatrix, KClusters, RegParam, ThresholdValue):
    #variables
    Centroids = np.empty([KClusters, DataMatrix.shape[1]], dtype=float) #corresponds to V in paper

    MembershipMatrix = np.empty([DataMatrix.shape[0], KClusters], dtype=float) #should be of shape i rows and k columns; corresponds to U in paper
    MatrixH = np.empty([DataMatrix.shape[0], KClusters], dtype=float) #derived from our s aux variable
    S = np.empty([DataMatrix.shape[0], KClusters], dtype=float) #holds our calulated s_ik
    TimeStep = 0

    print(type(MembershipMatrix))

    #initalization
    for rndx, row in enumerate(DataMatrix):
        for col in range(0, KClusters):
            MembershipMatrix[rndx][col] = float(1/KClusters) #this needs to be fixed to enforce constraint of U1 = 1 in cases of odd number of clusters like 3
            S[rndx][col] = 1

    Centroids = FindCentroids(DataMatrix, Centroids, S, MembershipMatrix)
