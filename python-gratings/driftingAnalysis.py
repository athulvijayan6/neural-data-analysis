#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-04 16:42:24
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-05 13:53:48

import numpy as np
import scipy.io
import osi

dataTargets = ['dataset/Mouse-A/', 'dataset/Mouse-B/', 'dataset/Mouse-C/', 'dataset/Mouse-D/', 'dataset/Mouse-E/']
for mouse in xrange(1):
    data = scipy.io.loadmat('../driftingGratings/'+ dataTargets[mouse] + 'Data.mat')
    data = data['Data']
    rawData = data[0, 0]['rawF']
    smoothData = data[0, 0]['dFF']
    stimuliSeq = data[0, 0]['StimSeq']
    cellData = np.zeros((smoothData.shape[0], stimuliSeq.size, 121))
    spikeRate = np.zeros((smoothData.shape[0], stimuliSeq.size, 2))
    OSI = np.zeros(smoothData.shape[0])
    DSI = np.zeros((smoothData.shape[0]))
    cirvar = np.zeros((smoothData.shape[0]))
    dircirvar = np.zeros((smoothData.shape[0]))
    for i in xrange(smoothData.shape[0]):
        for j in xrange(stimuliSeq.size):
            cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
        # sort the array
        cellData[i] = cellData[i, np.argsort(cellData[i, :, -1])]
        # Calculate spike rate response of each neuron as a real number
        spikeRate[i] = osi.calculateSpikeRate(cellData[i])
        OSI_n = DSI_n = np.zeros(10)
        for trial in xrange(10):
            trialData = spikeRate[i, trial:10:end, :]
            R_pref_ori = np.max(trialData[:, 1])
            theta_pref_ori = np.argmax(trialData[:, 1])
            # compute OSI
            theta_orth = ((theta_pref_ori + 3) % 16) + 1
            R_orth = trialData[theta_orth, 1]
            OSI_n[trial] = (R_pref_ori - R_orth)/(R_pref_ori + R_orth)
            # % Compute DSI
            theta_null = ((theta_pref_ori + 7) % 16) + 1
            R_null = trialData[theta_null, 1]
            DSI_n[trial] = (R_pref_ori - R_null)/(R_pref_ori + R_null)
        OSI[i] = np.mean(OSI_n)
        DSI[i] = np.mean(DSI_n)
        # % compute orientation variance
        cirvar[i] = osi.cirVar(spikeRate[i])
        # % compute directional variance
        dircirvar[i] = osi.dirCirVar(spikeRate[i])
