#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul Vijayan
# @Date:   2015-09-19 22:47:47
# @Last Modified by:   athul
# @Last Modified time: 2015-12-10 11:04:23
from __future__ import division
import numpy as np
import scipy.io
import osi.osi as osi
import curvefit.fit as fit
from scipy.cluster.vq import kmeans, whiten, vq
import matplotlib.pyplot as plt
from matplotlib import gridspec
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

    # ========================= Fit curve =============================
    numNeurons = spikeRate.shape[0]
    w0 = np.array([4, 50, np.pi/2, 1, 50, 3*np.pi/2, 1])
    what = np.zeros((numNeurons, w0.size))
    sse = np.zeros((numNeurons,))
    for n in xrange(numNeurons):
        what[n], sse[n] = fit.fitCurve(spikeRate[n, :, :], w0)


    # % ----------------------- k-means clustering ----------
    if True:
        # OSI, DSI, what[:, 3], what[:, 6]
        featureVec = np.vstack((np.abs(cirvar), np.abs(dircirvar), what[:, 1], what[:, 4], what[:, 3], what[:, 6]))
        featureVec = np.transpose(featureVec)
        whiteFeatureVec = whiten(featureVec)
        numClusters = 3
        centroids = kmeans(whiteFeatureVec, numClusters)[0]
        idx = vq(whiteFeatureVec, centroids)[0]

        # % Analyze the k means clustering results.
        globalMean = np.mean(whiteFeatureVec, 0)
        total_ss = 0
        for i in xrange(featureVec.shape[0]):
            total_ss = total_ss + np.linalg.norm(globalMean - whiteFeatureVec[i, :])
        between_ss = 0
        for i in xrange(numClusters):
            between_ss = between_ss + np.bincount(idx)[i]*np.linalg.norm(globalMean - centroids[i, :])
        clusterEfficiency = (between_ss/total_ss)*100
        print('The clustering is done with an efficiency of ' + str(clusterEfficiency))
        print('To see the definition of clustering efficiency, see the documentation')
    # -----------------------------------------------------------------------
    # % ----------------------- --------------------------------
    if True:
        # Plot k - means results on a 2D subspace of feature space
        oriData = {0:[], 1:[], 2:[]}
        dirData = {0:[], 1:[], 2:[]}
        names = {0:[], 1:[], 2:[]}
        j = 0
        for i in idx:
            oriData[i].append(featureVec[j, 0])
            dirData[i].append(featureVec[j, 1])
            names[i].append(str(j))
            j += 1
        fig, ax = plt.subplots()
        c = ['crimson', 'green', 'blue']
        for i in range(3):
            ax.scatter(oriData[i], dirData[i], label='class '+str(i), color = c[i], linewidth=8, alpha = 0.6)
        ax.set_title('Classification of simple and complex cells')
        ax.set_xlabel('Orientation index')
        ax.set_ylabel('Direction index')
    # ----------------------------------------------------------
    # ========================== Plot fits ============================
    if False:    
        for n in []:
            plt.figure()
            plt.subplot(111)
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])
            for trial in xrange(10):
                trialData = spikeRate[n, trial::10, :]
                ax0.plot(np.degrees(trialData[:, 1]), trialData[:, 0], linewidth=1, label='expe data of trial'+str(trial))
                res = trialData[:, 0] - fit.doubleGaussian(trialData[:, 1], what[n])
                ax1.plot(np.degrees(trialData[:, 1]), res, linewidth=1, label='trial '+str(trial))
            x = np.arange(0, 360, 5)
            y = fit.doubleGaussian(np.radians(x), what[n])
            ax1.plot(x, np.ones(x.shape), color='black')
            ax0.plot(x, y, linewidth=3.5, color='crimson', label='Fitted curve')

            ax0.set_title('Fitting of directional selectivity')
            ax0.set_xlabel(r'$\theta$')
            ax0.set_ylabel('spikerate')

            ax1.set_title('Residuals')
            ax1.set_xlabel(r'$\theta$')
            ax1.set_ylabel('residuals')
            plt.annotate('SSE: '+str(sse[n]), xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')
    
    plt.legend(loc='upper right')
    plt.show()