#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2015-09-04 15:23:53
# @Last Modified by:   Athul Vijayan
# @Last Modified time: 2015-09-06 17:54:47
import numpy as np

def calculateSpikeRate(conc, algo='average'): 
    '''the function computes the neuron response of a stimulus using the time-series. The function returns a real value which can be thought of as the spike rate during the stimulus period.'''
    conc = np.asarray(conc)
    spikeRate = np.zeros((conc.shape[0], 2))
    for i in xrange(conc.shape[0]):
        if algo == 'average':
            # This algorithm just computes the average Ca concentration
            # from the time of epoch to end of stimuli.The input expected is smoothed Ca data
            spikeRate[i, 0] = np.mean(conc[i, 80:])
            spikeRate[i, 1] = conc[i, -1]
    return spikeRate

def cirVar(spikeRate):
    '''Function calculates the L_ori from spike rate response of neuron to each angle and direction of stimuli'''
    cv = np.dot(spikeRate[:, 0],np.exp(2j*spikeRate[:, 1]))/np.sum(spikeRate[:, 0])
    return np.abs(cv)

def dirCirVar(spikeRate):
    '''Function calculates the L_dir from spike rate response of neuron to each angle and direction of stimuli'''
    dcv = np.dot(spikeRate[:, 0],np.exp(1j*spikeRate[:, 1]))/np.sum(spikeRate[:, 0])
    return np.abs(dcv)


if __name__=='__main__':
    print('Module for orientation sensitivity analysis of neurons in a drifting grating experiment')