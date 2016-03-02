#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-23 15:20:54
# @Last Modified by:   Athul
# @Last Modified time: 2016-03-01 14:33:26
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

def backtrack(c, s, d, a, X, Y, i, j):
    if (i == 0) or (j == 0):
        return np.array([])
    elif d[i, j] == 2:
        return [np.append(Z, X[i]) for Z in backtrack(c, s, d, a, X, Y, i-1, j-1)]
    elif d[i, j] == 1:
        return backtrack(c, s, d, a, X, Y, i-1, j)
    elif d[i, j] == 0:
        return backtrack(c, s, d, a, X, Y, i, j-1)

def lcss(X, Y, thres_sim=0.1, thres_rc=0.5, rho=1):
    m = Y.size
    n = X.size

    X = np.append(X, np.inf)
    Y = np.append(Y, np.inf)

    c = np.zeros((n+2, m+2)) # for storing running score
    s = np.zeros((n+2, m+2)) # for storing score
    d = np.zeros((n+2, m+2)) # For backtracking up = 0, left = 1, cross = 2
    a = np.zeros((n+2, m+2)) # for storing partial scores.
    p = min(m, n)
    for i in xrange(1, m+1):
        for j in xrange(1, n+1):
            xi = X[i]
            yj = Y[j]
            if sim(xi, yj) < thres_sim: #similar
                d[i, j] = 2
                c[i, j] = c[i-1, j-1] + (sim(xi, yj) - thres_sim)/(1 - thres_sim)
                s[i, j] = s[i-1, j-1]
            elif c[i-1, j] < c[i, j-1]:
                d[i, j] = 0
                c[i, j] = max(c[i-1, j] - rho, 0)
                if d[i-1, j] == 2:
                    s[i, j] = (a[i-1, j]*np.square(p) + np.square(c[i-1, j]))/np.square(p)
                else:
                    s[i, j] = s[i-1, j]
            else:
                d[i, j] = 1
                c[i, j] = max(c[i, j-1] - rho, 0)
                if d[i, j-1] == 2:
                    s[i, j] = (a[i, j-1]*np.square(p) + np.square(c[i, j-1]))/np.square(p)
                else:
                    s[i, j] = s[i, j-1]
            q = max(a[i-1, j-1], a[i-1, j], a[i, j-1])
            if (q == -1) and (d[i, j] == 2):
                a[i, j] = s[i-1, j-1]
            elif (c[i, j] < thres_rc):
                a[i, j] = -1
            else:
                a[i, j] = q
    return c, s, d, a




    




        
