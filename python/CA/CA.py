"""
Competitive Agglomeration Algorithm

Algorithm by:   Hichem Frigui and Raghu Krsihnapuram
Python Code:    Alain Garcia
"""

import numpy as np
import math

# Smallest number to be considered significant
EPSILON = 1.0e-6

# Largest number to be considered
INFINITY = 1.0e35

"""
Point(): Object used to store points
membership: A vector that holds membership for all clusters
distance:   A vector that holds distance from all clusters
dimension:  A vector that holds all the data for a certain point, separated by dimension
cluster:    Currently assigned cluster
"""
class Point():
    def __init__(self, membership, distance, dimension, cluster):
        self.membership = membership
        self.distance = distance
        self.dimension = dimension
        self.cluster = cluster

"""
The main Competitive Agglomeration function.
Takes in data, initialized centers, data details (vector_num, dimensions), and cluster number
Returns final number of cluster, final cluster center locations, and final point classifications
"""
def CA(in_data, in_centers,
        max_iterations=50, cluster_num=10, vector_num=60, dimensions=2,
        Eita_0=2.5, TAU=10, EXPINC=25, EXPDEC=35, 
        MAX_CENTER_DIFFERENCE=0.001, fuzzifier=2):
    
    # Initialize Data (from in_data)
    point_list = []
    for i in range(0, vector_num):
        point_list.insert(i,Point(np.zeros(cluster_num), 
                        np.zeros(cluster_num), np.zeros(dimensions), None))
    # Populate
    for i in range(0, dimensions):
        for j in range(0, vector_num):
            point_list[j].dimension[i] = in_data.iloc[j,i]

    # Initialize centerss (random data from in_centers)
    centers = in_centers.values

    # mimimum_points are used to determine the minimum number of points allowed in a cluster
    mimimum_points = np.zeros(max_iterations)
    mimimum_points2 = np.ones(max_iterations)
    mimimum_points2[int(max_iterations/2):max_iterations] *= 5

    # Keep track of old values 
    previous_cluster_num = cluster_num
    previous_centers = np.zeros((cluster_num, dimensions))
    for i in range(0, cluster_num):
        for j in range(0, dimensions):
            previous_centers[i][j] = centers[i][j]

    agglomeration_constant = 0
    iteration_num = 0

    # Run algorithm for maximum number of iterations
    while (iteration_num < max_iterations):

        # Calculate distance between each point and each cluster
        EuclideanDistance(point_list, centers, cluster_num, vector_num, dimensions, fuzzifier)
        
        # Assign points to closest cluster
        AssignPoints(point_list, cluster_num, vector_num)

        # If at least three iterations have not occurred, calculate membership using Fuzzy Algorithm
        if (iteration_num > 2):
            agglomeration_constant = AgglomerationConstant(point_list, cluster_num, iteration_num, vector_num, Eita_0, TAU, EXPINC, EXPDEC)
            Membership(point_list, cluster_num, agglomeration_constant, vector_num)
            cluster_num = ClusterNumber(point_list, cluster_num, iteration_num, vector_num, mimimum_points, mimimum_points2)
        else:
            FuzzyMembership(point_list, cluster_num, vector_num)
        
        # Use Fuzzy algorithms to calculate centers
        FuzzyCenters(point_list, centers, cluster_num, vector_num, dimensions)
        
        # Make sure the difference in center location has not changed too much
        centers_difference = 0
        for i in range(0, cluster_num):
            for j in range(0, dimensions):
                centers_difference += abs(previous_centers[i][j]-centers[i][j])
                
        if (centers_difference < MAX_CENTER_DIFFERENCE and previous_cluster_num == cluster_num):
            break
        
        previous_cluster_num = cluster_num

		# Copy the centers into the previous centers
        for i in range(0, cluster_num):
            for j in range(0, dimensions):
                previous_centers[i][j] = centers[i][j]
        
        iteration_num += 1 
    
    # Do a final calculation for distance, and a final point assignment
    EuclideanDistance(point_list, centers, cluster_num, vector_num, dimensions, fuzzifier)
    AssignPoints(point_list, cluster_num, vector_num)
    
    # Fill a list with classifications to return as output
    classification_list = []
    for i in range(0, vector_num):
        classification_list.append(point_list[i].cluster)

    # Fill a list with cluster centers to return as output
    center_list = []
    for i in range(0, cluster_num):
        center_list.append([])
        for j in range(0, dimensions):
            center_list[i].append(centers[i][j])

    return(cluster_num, center_list, classification_list)

# Calculate the distance from each point to each cluster
def EuclideanDistance(point_list, centers, cluster_num, vector_num, dimensions, fuzzifier):
    for i in range(0, vector_num):
        for j in range(0, cluster_num):
            point_list[i].distance[j] = 0
            for k in range(0, dimensions):
                temp = point_list[i].dimension[k]-centers[j][k]
                temp1 = 1
                for p in range(0, fuzzifier):
                    temp1 = temp1*temp
                point_list[i].distance[j] += abs(temp1)
            point_list[i].distance[j] = max(EPSILON, point_list[i].distance[j])

# Assign points to the closest avaiable cluster
def AssignPoints(point_list, cluster_num, vector_num):
    for i in range(0, vector_num):
        min = INFINITY
        for j in range(0, cluster_num):
            if point_list[i].distance[j] < min:
                min = point_list[i].distance[j]
                index = j
            point_list[i].cluster = index

# Calculates the agglomeration constant which controls the rate at which clusters are pruned
def AgglomerationConstant (point_list, cluster_num, iteration_num, vector_num, Eita_0, TAU, EXPINC, EXPDEC):
    objective_function, cardinality_sum = 0, 0
    cardinality = np.zeros(cluster_num)

    for i in range(0, cluster_num):
        cardinality[i] = 0
        for j in range(0, vector_num):
            temp = point_list[j].membership[i]
            cardinality[i] += temp
            objective_function += temp * temp * point_list[j].distance[i]
        cardinality_sum += cardinality[i] * cardinality[i]

    if (iteration_num < EXPINC):
        Exponent = (iteration_num - EXPINC)/TAU  
    elif (iteration_num < EXPDEC):
        Exponent = 0
    else:
        Exponent = (EXPDEC - iteration_num) / TAU
    
    alpha = Eita_0 * math.exp(Exponent) * objective_function/cardinality_sum
    
    return alpha

# Calculate the membership for all clusters for all points
# Membership contains a Fuzzy membership term and a membership bias term
def Membership(point_list, cluster_num, agglomeration_constant, vector_num):
    cardinality = np.zeros(cluster_num)
    
    # Calculate cardinality for each cluster (number of points in each cluster)
    for j in range(0, cluster_num):
        cardinality[j] = 0
        if (agglomeration_constant >= 0.0):
            for i in range(0, vector_num):
                cardinality[j] += point_list[i].membership[j]
    
    # Calculate inverse distance sum and average cardinality
    for i in range(0, vector_num):
        inverse_distance_sum, cardinality_average = 0.0, 0.0
        for j in range(0, cluster_num):
            inverse_distance_sum += 1.0/point_list[i].distance[j]
        for j in range(0, cluster_num):
            cardinality_average += cardinality[j]/point_list[i].distance[j]
        cardinality_average /= inverse_distance_sum

        for j in range(0, cluster_num):
            membership_fcm = 1.0 / (inverse_distance_sum * point_list[i].distance[j])
            membership_bias = (cardinality[j] - cardinality_average) / point_list[i].distance[j]
            point_list[i].membership[j] = membership_fcm + agglomeration_constant * membership_bias

            # Force membership to be between 0 and 1
            point_list[i].membership[j] = max(0, point_list[i].membership[j])
            point_list[i].membership[j] = min(1, point_list[i].membership[j])

        membership_sum = 0
        for j in range(0, cluster_num):
            membership_sum += point_list[i].membership[j]

        for j in range(0, cluster_num):
            point_list[i].membership /= membership_sum
        
        membership_sum = 0
        for j in range(0, cluster_num):
            membership_sum += point_list[i].membership[j]
        
        if membership_sum > 1.0e-10:
            for j in range(0, cluster_num):
                point_list[i].membership[j] /= membership_sum
        else:
            for j in range(0, cluster_num):
                point_list[i].membership[j] = 1 / cluster_num

# Calculate new number of clusters
def ClusterNumber(point_list, cluster_num, iteration_num, vector_num, mimimum_points, mimimum_points2):
    optimal_cluster_num = 0
    empty_clust = np.zeros(cluster_num)
    cardinality, cardinality2, cardinality3 = np.zeros(cluster_num), np.zeros(cluster_num), np.zeros(cluster_num)

    for i in range(0, cluster_num):
        empty_clust[i], cardinality[i], cardinality3[i] = 0, 0, 0
        for j in range(0, vector_num):
            cardinality[i] += point_list[j].membership[i]*point_list[j].membership[i]
            cardinality3[i] += point_list[j].membership[i]
            if point_list[j].cluster == i:
                cardinality2[i] += 1
        
        if (cardinality[i] <= mimimum_points[iteration_num]) or (cardinality2[i] <= mimimum_points2[iteration_num]):
            empty_clust[i] = 1

    for i in range(0, cluster_num):
        if empty_clust[i] == 0:
            for j in range(0, vector_num):
                point_list[j].membership[optimal_cluster_num] = point_list[j].membership[i]
                if point_list[j].cluster == i:
                    point_list[j].cluster = optimal_cluster_num
            optimal_cluster_num += 1

    return(optimal_cluster_num)

# Use Fuzzy algorithm to calculate membership
def FuzzyMembership(point_list, cluster_num, vector_num):
    for i in range(0, vector_num):
        inverse_distance_sum = 0
        for j in range(0, cluster_num):
            inverse_distance_sum +=  1 / point_list[i].distance[j]
        for j in range(0, cluster_num):
            temp = point_list[i].distance[j]
            point_list[i].membership[j] = 1 / (inverse_distance_sum*temp)

# Use Fuzzy algorithm to calculate centers
def FuzzyCenters(point_list, centers, cluster_num, vector_num, dimensions):
    cardinality = np.zeros(cluster_num)

    for j in range(0, cluster_num):
        cardinality[j] = 0
        for k in range(0, dimensions):
            centers[j][k] = 0

        for i in range(0, vector_num):
            temp = point_list[i].membership[j] * point_list[i].membership[j]

            cardinality[j] += temp
            for k in range(0, dimensions):
                centers[j][k] += temp * point_list[i].dimension[k]

        for k in range(0, dimensions):
            centers[j][k] /= cardinality[j] + EPSILON

# Initialize centers in an evenly spaced manner
# Currently not used because of the use of predefined pseudorandom center points
def InitializeCenters (centers, point_list, cluster_num, vector_num, dimensions):
    for k in range(0,dimensions):
        for j in range(0, cluster_num):
            index = j * vector_num/cluster_num
            centers[j][k] = point_list[index].dimension[k]