#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-21 11:41:30
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-23 13:08:47
from __future__ import division
import numpy as np
import scipy.io
import osi.osi as osi
import curvefit.fit as fit
import matplotlib.pyplot as plt
plt.style.use('ggplot')

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
    DSI = np.zeros(smoothData.shape[0])
    cirvar = np.zeros((smoothData.shape[0]))
    dircirvar = np.zeros((smoothData.shape[0]))
    for i in xrange(smoothData.shape[0]):
        for j in xrange(stimuliSeq.size):
            cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], stimuliSeq[j] )
        # sort the array
        cellData[i] = cellData[i, np.argsort(cellData[i, :, -1])]
        # Calculate spike rate response of each neuron as a real number
        spikeRate[i] = osi.calculateSpikeRate(cellData[i])

        cirvar[i], dircirvar[i], OSI[i], DSI[i] = np.abs(osi.computeAll(spikeRate[i]))

    # ======================= Correlation study ====================
    # Rearrange elements of spikeRate in decreasing circular variance

    spikeRate = spikeRate[np.argsort(cirvar)[::-1]]
    numNeurons = spikeRate.shape[0]
    corrMat = np.corrcoef(spikeRate[:, :, 0])

    # ********************** Plot correlation heat map *********
    fig, ax = plt.subplots()
    ax.pcolor(corrMat, cmap=plt.cm.Blues, alpha=0.8)
    # Format
    fig = plt.gcf()
    fig.set_size_inches(8, 11)
    # turn off the frame
    ax.set_frame_on(False)
    ax.grid(False)
    plt.title('Correlation heatmap')

    fig, ax = plt.subplots()
    corrThres = 0.6
    ax.pcolor(corrMat>corrThres, cmap=plt.cm.Blues, alpha=0.8)
    # Format
    fig = plt.gcf()
    fig.set_size_inches(8, 11)
    # turn off the frame
    ax.set_frame_on(False)
    ax.grid(False)
    plt.title('Correlation heatmap thresholded c > '+ str(corrThres))

    # **************** Plot response of selected cells *********

    selNeurons = np.array([1, 2, 3, 4, 5, 6, 7])
    numNeurons = selNeurons.size
    w0 = np.array([4, 50, np.pi/2, 1, 50, 3*np.pi/2, 1])
    what = np.zeros((numNeurons, w0.size))
    sse = np.zeros((numNeurons,))
    for n in xrange(numNeurons):
        what[n], sse[n] = fit.fitCurve(spikeRate[selNeurons[n], :, :], w0)
    fig, ax = plt.subplots()
    for n in xrange(numNeurons):
        x = np.arange(0, 360, 5)
        y = fit.doubleGaussian(np.radians(x), what[n])
        ax.plot(x, y, linewidth=3, label='Fitted curve of neuron '+ str(selNeurons[n]))
    ax.legend(loc='upper right')
    ax.set_title('Fitting of directional selectivity of multiple neurons')
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('spikerate')
    plt.show()