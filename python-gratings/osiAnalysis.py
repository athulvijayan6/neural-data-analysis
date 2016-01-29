#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-04 16:42:24
# @Last Modified by:   Athul
# @Last Modified time: 2015-11-17 11:49:42
from __future__ import division
import numpy as np
import scipy.io
import osi.osi as osi
from scipy.cluster.vq import kmeans, whiten, vq
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataTargets = ['dataset/Mouse-A/', 'dataset/Mouse-B/', 'dataset/Mouse-C/', 'dataset/Mouse-D/', 'dataset/Mouse-E/']
for mouse in xrange(1):
    data = scipy.io.loadmat('../driftingGratings/'+ dataTargets[mouse] + 'Data.mat')
    data = data['Data']
    rawData = data[0, 0]['rawF']
    smoothData = data[0, 0]['Spks']
    stimuliSeq = data[0, 0]['StimSeq']

    # To verify
    data = scipy.io.loadmat('../driftingGratings/'+ dataTargets[mouse] + '/Solutions/Ori.mat')
    data = data['Ori']
    OI_ref = data['OI'][0,0].reshape((65,))
    OSI_ref = data['OSI'][0,0].reshape((65,))
    
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

    # ==================== Verify results ===================
    if False:
        print(np.corrcoef(np.abs(cirvar), OSI_ref))   # Higher correlation required
        fig, ax = plt.subplots()
        x = np.array(xrange(OSI.size))
        ax.plot(x, OSI_ref, linewidth=2, label='reference')
        ax.plot(x, np.abs(cirvar), linewidth=2, label='estimated')

    # % ----------------------- k-means clustering ----------
    if False:
        featureVec = np.vstack((np.abs(cirvar), np.abs(dircirvar)))
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

    # ========================== Gaussian fit ===================================
    
    # % ====================== Plots =======================
    neuronId = 63

    # % ----------------------- --------------------------------
    if False:
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
        for i in range(3):
            ax.scatter(x = oriData[i], y = dirData[i], label='class '+str(i), alpha=0.6)
        ax.set_title('Classification of simple and complex cells')
        ax.set_xlabel('Orientation index')
        ax.set_title('Direction index')
    # ----------------------------------------------------------

    # % % 1. Plots each time-series data of length 120 with 
    # % % each trial of an experiment with seperate color
    if False:
        theta = 1
        fig, ax = plt.subplots()
        for j in xrange(10):
            x = [2*i/120 for i in xrange(121)]
            y = list(cellData[neuronId, 10*theta+j, :-2])
            ax.plot(x, y, label='trial '+str(trial))
        ax.set_title('Plots of time-series data with each trial in diffenrent color')
        ax.set_xlabel('time')
        ax.set_title('Spikes')
    # # % --------------------------------------------------------

    # % --------------------------3-----------------------------
    # % 3. Plots response vs angle of stimulus pattern.
    if False:
        fig, ax = plt.subplots()
        for trial in xrange(10):
            trialData = spikeRate[neuronId, trial::10, :]
            ax.plot(trialData[:, 1], trialData[:, 0], label='trial '+str(trial))
        ax.set_title('Neuron response (spike rate) Vs direction of stimuli')
        ax.set_xlabel('angle in radians')
        ax.set_title('Spike rate')
    # % --------------------------------------------------------


    # % --------------------------4-----------------------------
    # % 4. Polar plot of Directional selectivity
    # % Take a neuron, there will be 10 plots corresponding to each trial of the experiment
    # % plot each of the 10 plots in a figure to verify the similarity in response
    if True:
        plt.figure()
        ax = plt.subplot(111, polar=True)
        plt.grid(True)
        for trial in xrange(10):
            # % Normalize all trial response from 0 to 1
            normData = np.abs(spikeRate[neuronId, trial::10, 0])/np.max(np.abs(spikeRate[neuronId, trial::10, 0]));
            normData = np.append(normData, normData[0])
            try:
                avgData = avgData + normData
            except NameError:
                avgData = normData
            theta = spikeRate[neuronId, trial::10, 1]
            theta = np.append(theta, theta[0])
            ax.plot(theta, normData, linewidth=2, alpha=0.3)
        ax.plot(theta, avgData/10, linewidth=2, label='Average response', color='red')
        ax.set_title('Polar plot of Directional selectivity')
        ax.text(1.2, 0.8, 'DSI= %.2f' % np.abs(dircirvar[neuronId]))
        plt.legend(loc='upper right')
    # % --------------------------------------------------------

    # % --------------------------5-----------------------------
    # % 4. Polar plot of Orientation selectivity
    # % Take a neuron, there will be 10 plots corresponding to each trial of the experiment
    # % plot each of the 10 plots in a figure to verify the similarity in response

    # % There are observations for 16 angles, but there are only 8 orientations.
    # % We can treat them as two different observations. More data !
    if True:
        plt.figure()
        ax = plt.subplot(111, polar=True)
        plt.grid(True)
        del avgData
        for trial in xrange(10):
            # % Normalize all trial response from 0 to 1
            normData = np.abs(spikeRate[neuronId, trial::10, 0])/np.max(np.abs(spikeRate[neuronId, trial::10, 0]));
            try:
                avgData = avgData + normData
            except NameError:
                avgData = normData
            normData = np.append(normData, normData[0])
            theta = 2*spikeRate[neuronId, trial::10, 1]
            theta = np.append(theta, theta[0])
            ax.plot(theta, normData, linewidth=2, alpha=0.3)
        avgData = (avgData[:8] + avgData[8:])/2
        avgData = np.append(avgData, avgData[0])
        theta = np.append(theta[:8], theta[0])
        ax.plot(theta, avgData/10, linewidth=2, label='Average response', color='red')
        ax.set_xticklabels(['0', '', '45', '', '90', '', '135', ''])
        ax.set_title('Polar plot of Orientation selectivity')
        ax.text(1, 0.8, 'neuron with OSI= %.2f' % np.abs(cirvar[neuronId])+'\nPreffered orientation= %.2f' % np.angle(cirvar[neuronId], deg=True)+' deg')

    plt.legend(loc='upper right')
    plt.show()    
    # % ********************** Plots ************************

