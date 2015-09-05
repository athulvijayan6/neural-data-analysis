#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2015-09-04 15:23:53
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-05 13:48:47
import numpy as np

def calculateSpikeRate(conc, algo='average'):
    print(type(conc))
    conc = np.asarray(conc)
    spikeRate = np.zeros((conc.shape[0], 2))
    for i in xrange(conc.shape[0]):
        if algo == 'average':
            # This algorithm just computes the average Ca concentration
            # from the time of epoch to end of stimuli.The input expected is smoothed Ca data
            spikeRate[i, 1] = np.mean(conc[i, 80:])
            spikeRate[i, 2] = conc[i, -1]
            return spikeRate

def cirVar(spikeRate):
    spikeRate = np.asarray(spikeRate)
    cv = np.dot(spikeRate[:, 1],np.exp(2j*spikeRate[:, 2]))/np.sum(spikeRate[:, 1])
    return cv

def dirCirVar(spikeRate):
    spikeRate = np.asarray(spikeRates)
    dcv = np.dot(spikeRate[:, 1],np.exp(1j*spikeRate[:, 2]))/np.sum(spikeRate[:, 1])
    return dcv


if __name__=='__main__':
    print('Module for orientation sensitivity analysis of neurons in a drifting grating experiment')