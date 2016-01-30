#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-21 11:41:30
# @Last Modified by:   athul
# @Last Modified time: 2015-12-10 11:04:30
from __future__ import division
import numpy as np
import scipy.io
import osi.osi as osi
import curvefit.fit as fit
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataTargets = ['../datasets/driftingGratings/Mouse-A/', '../datasets/driftingGratings/Mouse-B/', '../datasets/driftingGratings/Mouse-C/', '../datasets/driftingGratings/Mouse-D/', '../datasets/driftingGratings/Mouse-E/']
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
    spikeRate = spikeRate[np.argsort(np.abs(cirvar))[::-1]]

    # ======================= Correlation study ====================
    # Rearrange elements of spikeRate in decreasing circular variance
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
    corrThres = 0.65
    ax.pcolor(corrMat>corrThres, cmap=plt.cm.Blues, alpha=0.8)
    # Format
    fig = plt.gcf()
    fig.set_size_inches(8, 11)
    # turn off the frame
    ax.set_frame_on(False)
    ax.grid(False)
    plt.title('Correlation heatmap thresholded c > '+ str(corrThres))

    # **************** Plot response of selected cells *********
    # [7, 10, 16, 23, 28]
    # [10, 16,23, 26, 28]
    # [14, 21, 30, 40, 44, 49, 56]
    # [16, 23, 26, 28]
    # [19, 27, 30, 33, 37, 38, 44, 47, 49, 50, 56, 58]
    # [21, 25, 30, 33, 36, 37, 38, 40, 44, 49, 50, 56, 58, 63]
    # [30, 33, 36, 37, 38, 40, 41, 44, 49, 50, 56, 58, 60, 63]
    selNeurons = np.array([30, 33, 36, 37, 38, 40, 41, 44, 49, 50, 56, 58, 60, 63])
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
    ax.set_title('Response curves of neurons showed high correlation (>0.7) ')
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('spikerate')
    plt.text(1,1, 'Neurons '+str(selNeurons)+' are tuned similarly, as their preferred orientation is same',)
    plt.show()