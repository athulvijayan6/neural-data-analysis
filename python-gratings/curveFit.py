#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul Vijayan
# @Date:   2015-09-19 22:47:47
# @Last Modified by:   Athul
# @Last Modified time: 2016-05-02 12:52:57
from __future__ import division
import numpy as np
import scipy.io
import datetime
import osi.osi as osi
import curvefit.fit as fit
from scipy.cluster.vq import kmeans, whiten, vq
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

plotDir = '../plots/'

dataTargets = ['../datasets/driftingGratings/Mouse-A/', '../datasets/driftingGratings/Mouse-B/', '../datasets/driftingGratings/Mouse-C/', '../datasets/driftingGratings/Mouse-D/', '../datasets/driftingGratings/Mouse-E/']
for mouse in xrange(1):
    data = scipy.io.loadmat(dataTargets[mouse] + 'Data.mat')
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
    spikeRate = spikeRate[np.argsort(np.abs(cirvar))[::-1]]

    # # ========================= Fit orientation curve =============================
    numNeurons = spikeRate.shape[0]
    w0_ori = np.array([4, 50, np.pi/2, 1])
    what_ori = np.zeros((numNeurons, w0_ori.size))
    sse_ori = np.zeros((numNeurons,))
    for n in xrange(numNeurons):
        spikeRate_ori = spikeRate
        spikeRate_ori[n, 80:, 1] = spikeRate[n, :80, 1]
        what_ori[n], sse_ori[n] = fit.fitCurve(spikeRate_ori[n, :80, :], w0_ori, func=fit.uniGaussian)

    # ========================= Fit direction curve =============================
    numNeurons = spikeRate.shape[0]
    w0_dir = np.array([4, 50, np.pi/2, 1, 50, 3*np.pi/2, 1])
    what_dir = np.zeros((numNeurons, w0_dir.size))
    sse_dir = np.zeros((numNeurons,))
    for n in xrange(numNeurons):
        what_dir[n], sse_dir[n] = fit.fitCurve(spikeRate[n, :, :], w0_dir, func=fit.doubleGaussian)

    # ========================== Plot fits ============================
    # Fit of orientation responses
    if True:    
        for n in [3]:
            fig = plt.figure()
            plt.subplot(111)
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])
            for trial in xrange(10):
                trialData = spikeRate[n, trial:80:10, :]
                ax0.plot(np.degrees(trialData[:, 1]), trialData[:, 0], linewidth=1)
                res = trialData[:, 0] - fit.uniGaussian(trialData[:, 1], what_ori[n])
                ax1.plot(np.degrees(trialData[:, 1]), res, linewidth=1)
            x = np.arange(0, 180, 5)
            y = fit.uniGaussian(np.radians(x), what_ori[n])
            ax1.plot(x, np.ones(x.shape), color='black')
            ax0.plot(x, y, linewidth=3, color='crimson', label='Fitted curve')

            ax0.set_title('Fitting of orientation selectivity')
            ax0.set_xlabel(r'$\theta$')
            ax0.set_ylabel('spikerate')

            ax1.set_title('Residuals')
            ax1.set_xlabel(r'$\theta$')
            ax1.set_ylabel('residuals')
            plt.annotate('SSE: '+str(sse_ori[n]), xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')
            now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            fig.savefig(plotDir+'gratings_crvefit_osi_n_'+str(n)+'_'+now+'.pdf')

    # Fit of direction responses
    if False:    
        for n in [0]:
            fig = plt.figure()
            plt.subplot(111)
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])
            for trial in xrange(10):
                trialData = spikeRate[n, trial::10, :]
                ax0.plot(np.degrees(trialData[:, 1]), trialData[:, 0], linewidth=1)
                res = trialData[:, 0] - fit.doubleGaussian(trialData[:, 1], what_dir[n])
                ax1.plot(np.degrees(trialData[:, 1]), res, linewidth=1)
            x = np.arange(0, 360, 5)
            y = fit.doubleGaussian(np.radians(x), what_dir[n])
            ax1.plot(x, np.ones(x.shape), color='black')
            ax0.plot(x, y, linewidth=3.5, color='crimson', label='Fitted curve')

            ax0.set_title('Fitting of directional selectivity')
            ax0.set_xlabel(r'$\theta$')
            ax0.set_ylabel('spikerate')

            ax1.set_title('Residuals')
            ax1.set_xlabel(r'$\theta$')
            ax1.set_ylabel('residuals')
            plt.annotate('SSE: '+str(sse_dir[n]), xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')
            now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            fig.savefig(plotDir+'gratings_crvefit_n_'+str(n)+'_'+now+'.pdf')
    
    plt.legend(loc='upper right')
    plt.show()