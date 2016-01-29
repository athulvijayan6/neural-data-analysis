#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul Vijayan
# @Date:   2015-09-19 22:47:47
# @Last Modified by:   Athul
# @Last Modified time: 2015-11-17 11:53:18
from __future__ import division
import numpy as np
import scipy.io
import osi.osi as osi
import curvefit.fit as fit
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

dataTargets = ['dataset/Mouse-A/', 'dataset/Mouse-B/', 'dataset/Mouse-C/', 'dataset/Mouse-D/', 'dataset/Mouse-E/']
for mouse in xrange(1):
    data = scipy.io.loadmat('../driftingGratings/'+ dataTargets[mouse] + 'Data.mat')
    data = data['Data']
    rawData = data[0, 0]['rawF']
    smoothData = data[0, 0]['Spks']
    stimuliSeq = data[0, 0]['StimSeq']
    cellData = np.zeros((smoothData.shape[0], stimuliSeq.size, 121))
    spikeRate = np.zeros((smoothData.shape[0], stimuliSeq.size, 2))
    OSI = np.zeros(smoothData.shape[0])
    DSI = np.zeros(smoothData.shape[0])
    cirvar = np.zeros((smoothData.shape[0]), dtype=np.complex64)
    dircirvar = np.zeros((smoothData.shape[0]), dtype=np.complex64)
    for i in xrange(smoothData.shape[0]):
        for j in xrange(stimuliSeq.size):
            cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
        # sort the array
        cellData[i] = cellData[i, np.argsort(cellData[i, :, -1])]
        # Calculate spike rate response of each neuron as a real number
        spikeRate[i] = osi.calculateSpikeRate(cellData[i])
        cirvar[i], dircirvar[i], OSI[i], DSI[i] = osi.computeAll(spikeRate[i])

    # ========================= Fit curve =============================
    numNeurons = spikeRate.shape[0]
    w0 = np.array([4, 50, np.pi/2, 1, 50, 3*np.pi/2, 1])
    what = np.zeros((numNeurons, w0.size))
    sse = np.zeros((numNeurons,))
    for n in xrange(numNeurons):
        what[n], sse[n] = fit.fitCurve(spikeRate[n, :, :], w0)

    # ========================== Plot fits ============================
    
    for n in [49]:
        plt.figure()
        plt.subplot(111)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])
        for trial in xrange(10):
            trialData = spikeRate[n, trial::10, :]
            ax0.plot(np.degrees(trialData[:, 1]), trialData[:, 0], alpha=0.5, linewidth=1, label='expe data of trial'+str(trial))
            res = trialData[:, 0] - fit.doubleGaussian(trialData[:, 1], what[n])
            ax1.plot(np.degrees(trialData[:, 1]), res, linewidth=1, label='trial '+str(trial))
        x = np.arange(0, 360, 5)
        y = fit.doubleGaussian(np.radians(x), what[n])
        ax1.plot(x, np.ones(x.shape), color='black')
        ax0.plot(x, y, linewidth=3.5, color='crimson', label='Fitted curve')

        ax0.set_title('Fitting of directional selectivity using double Gaussian')
        ax0.set_xlabel(r'$\theta$')
        ax0.set_ylabel('spikerate')

        ax1.set_title('Residuals')
        ax1.set_xlabel(r'$\theta$')
        ax1.set_ylabel('residuals')
        plt.annotate('SSE: '+str(sse[n]), xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')

        plt.legend(loc='upper right')
        plt.show()