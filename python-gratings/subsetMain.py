#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-04 14:55:40
# @Last Modified by:   athul
# @Last Modified time: 2016-02-05 00:05:28
from __future__ import division
import numpy as np
import scipy.io
import osi.osi as osi
import curvefit.fit as fit
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataTargets = ['../datasets/driftingGratings/Mouse-A/', '../datasets/driftingGratings/Mouse-B/', '../datasets/driftingGratings/Mouse-C/', '../datasets/driftingGratings/Mouse-D/', '../datasets/driftingGratings/Mouse-E/']
for mouse in xrange(1):
    data = scipy.io.loadmat(dataTargets[mouse] + 'Data.mat')
    data = data['Data']
    rawData = data[0, 0]['rawF']
    smoothData = data[0, 0]['Spks']
    stimuliSeq = data[0, 0]['StimSeq']
    stimuli = np.unique(stimuliSeq)
    numNeurons = smoothData.shape[0]
    avgData = np.zeros((stimuli.size, 120, numNeurons))
    cellData = np.zeros((numNeurons, stimuliSeq.size, 121))
    spikeRate = np.zeros((numNeurons, stimuliSeq.size, 2))

    cirvar = np.zeros((numNeurons), dtype=np.complex64)
    dircirvar = np.zeros((numNeurons), dtype=np.complex64)

    for i in xrange(numNeurons):
        for j in xrange(stimuliSeq.size):
            theta = stimuliSeq[j]
            cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
            avgData[stimuli == theta, :, i] += smoothData[i, 120*j: 120*(j+1)]
        avgData = avgData/stimuli.size

        # Calculate spike rate response of each neuron as a real number
        spikeRate[i] = osi.calculateSpikeRate(cellData[i])
        cirvar[i], dircirvar[i], _, _ = osi.computeAll(spikeRate[i])

    if False:
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.2, 0.8, 0.6])
        for n in xrange(numNeurons):
            ax.plot(xrange(120), avgData[0, :, n], label='n = '+str(n))
        ax.set_title('Average response of all neurons to a particular orientation')
        ax.set_xlabel('time')
        ax.set_ylabel('Spike rate')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=12)
    for i in xrange(stimuli.size):
        try:
            allResponses = np.vstack((allResponses, avgData[i]))
        except:
            allResponses = avgData[i]
    allResponses = preprocessing.scale(allResponses)
    pca = PCA()
    pca.fit(allResponses)
    fig, ax = plt.subplots()
    varDist = pca.explained_variance_ratio_
    ax.plot(xrange(len(varDist)), varDist)
    plt.show()











