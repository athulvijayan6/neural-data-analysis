#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-10 12:40:50
# @Last Modified by:   Athul
# @Last Modified time: 2016-05-02 12:23:23
import numpy as np
from scipy.optimize import leastsq

def gaussianCurve(x, a, b, c):
    '''Implements gaussian curve function
    f(x) = a exp(-(x-b)^2/2c^2)
    '''
    inner = np.square(x-b)/(2*np.square(c))
    return a*np.exp(-inner)

def gaussianPrime(Rtrue, theta, params):
    c, Rp, theta_1, sigma = list(params)
    dc = da = db = dsigma = 0
    for i in xrange(Rtrue.size):
        gaus = gaussianCurve(theta[i], Rp, theta_1, sigma)
        dc += 2*(Rtrue[i] - c - gaus)
        da += 2*(Rtrue[i] - c - gaus)*(-gaus)/Rp
        db += 2*(Rtrue[i] - c - gaus)*(-gaus)*(theta[i] - theta_1)/np.square(sigma)
        dsigma += 2*(Rtrue[i] - c - gaus)*(-gaus)*np.square(theta[i] - theta_1)/np.power(sigma, 3)
    return np.array([dc, da, db, dsigma])

def doubleGaussian(x, params):
    return params[0] + gaussianCurve(x, params[1], params[2], params[3]) + gaussianCurve(x, params[4], params[5], params[6])

def uniGaussian(x, params):
    return params[0] + gaussianCurve(x, params[1], params[2], params[3])

def residuals(params, R, theta, func=doubleGaussian):
    Rhat = func(theta, params);
    return (R - Rhat)

def fitCurve(spikeRate, w0, func=doubleGaussian):
    theta = spikeRate[:, 1]
    R = spikeRate[:, 0]
    what = leastsq(residuals, w0, args=(R, theta, func))[0]
    sse = sum(np.square(residuals(what, R, theta, func)))
    return what, sse


if __name__=='__main__':
    print('Module for general curve fitting')