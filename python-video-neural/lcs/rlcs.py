# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-03-07 18:04:47
# @Last Modified by:   Athul
# @Last Modified time: 2016-03-08 12:25:22
import numpy as np
from math import floor

def distance(x, y, measure="euclidean"): # distance measure.
    if (measure == "euclidean"):
        return np.linalg.norm(x - y)**2
    elif (measure == "cubic"):
        return np.power((x - y), 3)  # TODO - doesn't work for multi dim
    else:
        print("Invalid measure! ")

def backtrack(X, Y, score, diag, cost):
    # Find max and min of score matrix
    maxScore = -1*np.inf
    for i in xrange(score.shape[0]):
        for j in xrange(score.shape[1]):
            if (maxScore < score[i, j]):
                maxScore = score[i, j]
                x, y = i, j         # indices of max score
    segLen = 0
    segment = np.array([[0, 0, 0]])   # stores ref position, query position, color for each point
    while(x != 0) and (y != 0):
        i = segLen
        segment[i][0] = x - 1
        segment[i][1] = y - 1
        segment[i][2] = cost[x, y]
        segment = np.vstack((segment, np.array([0, 0, 0])))  # make row for next point
        segLen += 1
        # go to new point        
        if (diag[x, y] == 1):
            x -= 1
            y -= 1
        elif (diag[x, y] == 2):
            x -=1
        else:
            y -= 1
    segment = segment[:-1]          # One row added but while loop exited. remove that
    segment = np.flipud(segment)    # We were backtracking.
    return segment
    

def rlcs(X, Y, tau_dist=0.005, delta=0.5):
    # if the input is vector (1D problem), count the elements
    # if not, count number of rows. Each point is a row and each column is a dim
    m = X.shape[0]   # Expect matrix/ multidimensional input
    n = Y.shape[0]
    if (m == 1):     # If input is a vector, take it as 1D array
        m = X.size
    if (n == 1):
        n = Y.size

    # find min distance and max distance
    maxDist = minDist = distance(X[0], Y[0])
    for i in xrange(m):
        for j in xrange(n):
            dist = distance(X[i], Y[j])
            if dist > maxDist:
                maxDist = dist
            if dist < minDist:
                minDist = dist
    # Initialize matrices for dynamic programming
    cost = np.zeros((m+2, n+2)) # for storing running score
    score = np.zeros((m+2, n+2)) # for storing score
    diag = np.zeros((m+2, n+2)) # For backtracking cross = 0, up = 2, left = 3
    partial = np.zeros((m+2, n+2)) # for storing partial scores.
    p = min(m, n)

    loop_count = 0
    for i in xrange(1, m+1):
        for j in xrange(1, n+1):
            xi = X[i-1]
            yj = Y[j-1]
            # find the distance
            dist = (distance(xi, yj) - minDist)/(maxDist - minDist)
            # take dp action
            if dist < tau_dist:
                diag[i, j] = 1 # cross arrow
                cost[i, j] = cost[i-1, j-1] + (1 - dist/tau_dist)
                score[i, j] = score[i-1, j-1]
            elif (cost[i-1, j] > cost[i, j-1]):
                diag[i, j] = 2 # up arrow
                cost[i, j] = cost[i-1, j] - delta
                if cost[i, j] < 0:
                    cost[i, j] = 0
                if (diag[i-1, j] == 1):
                    score[i, j] = (partial[i-1, j]*np.square(p) + np.square(cost[i-1, j]))/np.square(p)
                else:
                    score[i, j] = score[i-1, j]
            else:
                diag[i, j] = 3 # left arrow
                cost[i, j] = cost[i, j-1] - delta
                if cost[i, j] < 0:
                    cost[i, j] = 0
                if (diag[i, j-1] == 1):
                    score[i, j] = (partial[i, j-1]*np.square(p) + np.square(cost[i-1, j]))/np.square(p)
                else:
                    score[i, j] = score[i, j-1]

            a = partial[i-1, j-1]
            b = partial[i-1, j]
            c = partial[i, j-1]
            mabc = max(a, b, c)

            if (mabc == -1) and (diag[i, j] == 1):
                partial[i, j] = score[i-1, j-1]
            # TODO Confusion here. Is it 0.5 or delta?
            elif (cost[i, j] <= delta) and (diag[i, j] != 1):
                partial[i, j] = -1
            else:
                partial[i, j] = mabc
            loop_count += 1

    lri = score.shape[0] - 1; # last row index of score matrix
    lci = score.shape[1] - 1; # last column index of score matrix

    for i in xrange (lri):
        # last row
        if (diag[lri-1, i] == 1):
            score[i, lci] = (partial[i, lci-1]*np.square(p) + np.square(cost[i, lci-1]))/np.square(p)
        else:
            score[i, lci] = score[i, lci-1]
        diag[lri, i] = 3                    # left

    for i in xrange (lci):
        # last row
        if (diag[lri-1, i] == 1):
            score[lri, i] = (partial[lri-1, i]*np.square(p) + np.square(cost[lri-1, i]))/np.square(p)
        else:
            score[lri, i] = score[lri-1, i]
        diag[lri, i] = 2                    # up

    score[lri, lci] = score[lri-1, lci]
    diag[lri, lci] = 2

    return score, diag, cost

if __name__ == '__main__':
    print(''' This module is not intended to run from interpreter.
        Instead, call the functions from your main script.
        from lcs import rlcs as rlcs

        score, diag, cost = rlcs.rlcs(X, Y, tau_dist,  delta)
        segment = rlcs.backtrack(X, Y, score, diag, cost)''')