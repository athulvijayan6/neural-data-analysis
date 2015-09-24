#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2015-09-04 15:23:53
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-24 11:13:59
import numpy as np

def calculateSpikeRate(conc, algo='average'): 
    '''the function computes the neuron response of a stimulus using the time-series. The function returns a real value which can be thought of as the spike rate during the stimulus period.'''
    conc = np.asarray(conc)
    spikeRate = np.zeros((conc.shape[0], 2))
    for i in xrange(conc.shape[0]):
        if algo == 'average':
            # This algorithm just computes the average Ca concentration
            # from the time of epoch to end of stimuli.The input expected is smoothed Ca data
            spikeRate[i, 0] = np.mean(conc[i, 80:-1])
            spikeRate[i, 1] = conc[i, -1]
    return spikeRate

def cirVar(spikeRate, angleUnit='radians'):
    '''Function calculates the L_ori from spike rate response of neuron to each angle and direction of stimuli'''
    if angleUnit != 'radians':
        spikeRate[:, -1] = np.radians(spikeRate[:, -1])
    cv = np.dot(spikeRate[:, 0],np.exp(2j*spikeRate[:, 1]))/np.sum(spikeRate[:, 0])
    return cv

def dirCirVar(spikeRate, angleUnit='radians'):
    '''Function calculates the L_dir from spike rate response of neuron to each angle and direction of stimuli'''
    if angleUnit != 'radians':
        spikeRate[:, -1] = np.radians(spikeRate[:, -1])
    dcv = np.dot(spikeRate[:, 0],np.exp(1j*spikeRate[:, 1]))/np.sum(spikeRate[:, 0])
    return dcv

def computeAll(spikeRate):
    OSI_n = DSI_n = np.zeros((10))
    for trial in xrange(10):
        trialData = spikeRate[trial::10, :]
        R_pref_ori = np.max(trialData[:, 0])
        theta_pref_ori = np.argmax(trialData[:, 0])
        # compute OSI
        theta_orth = (theta_pref_ori + 4) % 16
        R_orth = trialData[theta_orth, 0]
        OSI_n[trial] = (R_pref_ori - R_orth)/(R_pref_ori + R_orth)
        # % Compute DSI
        theta_null = (theta_pref_ori + 8) % 16
        R_null = trialData[theta_null, 0]
        DSI_n[trial] = (R_pref_ori - R_null)/(R_pref_ori + R_null)
    OSI= np.mean(OSI_n)
    DSI = np.mean(DSI_n)
    # % compute orientation variance
    cirvar = cirVar(spikeRate)
    # % compute directional variance
    dircirvar = dirCirVar(spikeRate)
    return cirvar, dircirvar, OSI, DSI

if __name__=='__main__':
    print('Module for orientation sensitivity analysis of neurons in a drifting grating experiment')