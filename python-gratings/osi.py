#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2015-09-04 15:23:53
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-05 16:56:38
import numpy as np

def calculateSpikeRate(conc):
    algorithm = 1
    conc = np.asarray(conc)
    spikeRate = np.zeros((conc.shape[0], 2))
    for i in xrange(conc.shape[0]):
        if algorithm == 1:
            # This algorithm just computes the average Ca concentration
            # from the time of epoch to end of stimuli.The input expected is smoothed Ca data
            spikeRate[i, 0] = np.mean(conc[i, 80:])
            spikeRate[i, 1] = conc[i, -1]
    return spikeRate

def cirVar(spikeRate):
    cv = np.dot(spikeRate[:, 0],np.exp(2j*spikeRate[:, 1]))/np.sum(spikeRate[:, 0])
    return np.abs(cv)

def dirCirVar(spikeRate):
    dcv = np.dot(spikeRate[:, 0],np.exp(1j*spikeRate[:, 1]))/np.sum(spikeRate[:, 0])
    return np.abs(dcv)


if __name__=='__main__':
    print('Module for orientation sensitivity analysis of neurons in a drifting grating experiment')