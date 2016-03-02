# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-03-01 14:29:54
# @Last Modified by:   Athul
# @Last Modified time: 2016-03-01 14:59:33
import numpy as np
def sim(x, y, measure="euclidean"):
    if (x == np.inf) or (y == np.inf):
        return 0
    if measure == "euclidean":
        return np.square(x - y)
    if measure == "cubic":
        return 1 - np.power((x - y), 3)
    else:
        print("Invalid measure! ")

def backtrack(X, Y, C, i, j):
    if (i == 0) or (j == 0):
        return np.array([])
    elif (sim(xi, yj) < tau_sim):


def lcss(X, Y, tau_sim=0.1, tau_rc=0.5, rho=1):
    n = X.size
    m = Y.size

    X = np.insert(X, np.inf, 0)
    Y = np.insert(Y, np.inf, 0)

    c = np.zeros((n+1, m+1))        # for storing running score
    s = np.zeros((n+1, m+1))        # for storing score
    d = np.zeros((n+1, m+1))        # For backtracking up = 0, left = 1, cross = 2
    a = np.zeros((n+1, m+1))        # for storing partial scores.
    p = min(m, n)
    for i in xrange(1, m+1):
        for j in xrange(1, n+1):
            xi = X[i]
            yj = Y[j]
            if sim(xi, yj) < tau_sim:       #similar
                d[i, j] = 2
                c[i, j] = c[i-1, j-1] + 1
            elif c[i-1, j] < c[i, j-1]:     # Left cell has more length
                d[i, j] = 1
                c[i, j] = c[i, j-1]
            else:                           # Right cell has more length
                d[i, j] = 0
                c[i, j] = c[i-1, j]
    return c, s, d, a