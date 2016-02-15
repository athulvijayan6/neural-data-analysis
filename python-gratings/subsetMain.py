#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-04 14:55:40
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-15 15:02:13
from __future__ import division
import datetime
import numpy as np
import scipy.io
import osi.osi as osi
import curvefit.fit as fit
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

dataTargets = ['../datasets/driftingGratings/Mouse-A/', '../datasets/driftingGratings/Mouse-B/', '../datasets/driftingGratings/Mouse-C/', '../datasets/driftingGratings/Mouse-D/', '../datasets/driftingGratings/Mouse-E/']
plotDir = '../plots/'
for mouse in xrange(1):
    data = scipy.io.loadmat(dataTargets[mouse] + 'Data.mat')
    data = data['Data']
    rawData = data[0, 0]['rawF']
    smoothData = data[0, 0]['Spks']
    stimuliSeq = data[0, 0]['StimSeq']
    stimuli = np.unique(stimuliSeq)
    numNeurons = smoothData.shape[0]
    # Do PCA to find eigen neurons
    numComps = [16]
    err = np.array([])
    # 
    for numEigNeurons in numComps:
        # ************************* Do PCA *********************************
        allData = np.transpose(smoothData)
        pca = PCA(n_components=numEigNeurons)
        pca.fit(allData)
        # ************************* PCA ends *********************************
        # Plot variance degradation
        varDist = pca.explained_variance_ratio_
        if False:
            fig, ax = plt.subplots()
            ax.plot(xrange(len(varDist)), varDist)
            ax.set_xlabel('Number of components')
            ax.set_ylabel('variance explained')
            ax.set_title('Principal component analysis')
            now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            fig.savefig(plotDir+'subsetbMain_pcaPlot_'+now+'.pdf')
        # compute reconstruction errors
        # Transform to new dimensions
        eigRes = pca.transform(allData)

        # Reconstruct back to original coordinates to find error
        allDataHat = pca.inverse_transform(eigRes)
        e = np.linalg.norm(allData - allDataHat)
        err = np.append(err, e)

        # **************** Do OSI analysis on Eigen neurons *******************

        eigRes = np.transpose(eigRes)
        # eigRes contain responses of 'Eigen neurons' - Neuron response in eigen space
        eigCellData = np.zeros((numEigNeurons, stimuliSeq.size, 121))
        eigSpikeRate = np.zeros((numEigNeurons, stimuliSeq.size, 2))

        eigCirvar = np.zeros((numEigNeurons), dtype=np.complex64)
        eigDircirvar = np.zeros((numEigNeurons), dtype=np.complex64)
        for i in xrange(numEigNeurons):
            for j in xrange(stimuliSeq.size):
                theta = stimuliSeq[j]
                eigCellData[i, j] = np.append(eigRes[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
            # Calculate spike rate response of each neuron as a real number
            eigSpikeRate[i] = osi.calculateSpikeRate(eigCellData[i])
            eigSpikeRate[i] = preprocessing.scale(eigSpikeRate[i])
            eigCirvar[i], eigDircirvar[i], _, _ = osi.computeAll(eigSpikeRate[i])
        # ************************ OSI analysis ends ****************************
        
        # scatter plot preferred orientation of each eigen neurons
        if False:
            x = np.angle(eigCirvar)
            y = numEigNeurons*np.ones(x.shape)
            if 'figPrefOri' not in locals():
                figPrefOri, axPrefOri = plt.subplots()
                axPrefOri.set_xlabel('Orientation')
                axPrefOri.set_ylabel('Number of principal components')
                axPrefOri.set_title('Orientation of eigen neurons for different components')
            axPrefOri.scatter(x, y, s=30, color='red')

            if numEigNeurons == numComps[-1]:
                now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
                figPrefOri.savefig(plotDir+'subsetbMain_prefOri_'+now+'.pdf')

        # Fit double gaussian for each eigen neuron
        w0 = np.array([4, 50, np.pi/2, 1, 50, 3*np.pi/2, 1])
        what = np.zeros((numEigNeurons, w0.size))
        sse = np.zeros((numEigNeurons,))
        for n in xrange(numEigNeurons):
            what[n], sse[n] = fit.fitCurve(eigSpikeRate[n, :, :], w0)
        # Plot double gaussians for each eigen neuron
        if False:
            fig, ax = plt.subplots()
            for n in xrange(numEigNeurons):
                x = np.arange(0, 360, 5)
                y = fit.doubleGaussian(np.radians(x), what[n])
                ax.plot(x, y, linewidth=2, label='Neuron '+str(n))
            ax.set_title('Fitting of directional selectivity for eigen neurons')
            ax.set_xlabel(r'$\theta$')
            ax.set_ylabel('spikerate')
            now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            fig.savefig(plotDir+'subsetbMain_dgauss_'+now+'.pdf')

        if True:
            n = 0
            fig = plt.figure()
            plt.subplot(111)
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])
            for trial in xrange(10):
                trialData = eigSpikeRate[n, trial::10, :]
                ax0.plot(np.degrees(trialData[:, 1]), trialData[:, 0], linewidth=1, label='trial'+str(trial))
                res = trialData[:, 0] - fit.doubleGaussian(trialData[:, 1], what[n])
                ax1.plot(np.degrees(trialData[:, 1]), res, linewidth=1)
            # x = np.arange(0, 360, 5)
            # y = fit.doubleGaussian(np.radians(x), what[n])
            # ax1.plot(x, np.zeros(x.shape), color='black')
            # ax0.plot(x, y, linewidth=3.5, color='crimson', label='Fitted curve')

            ax0.set_title('Fitting of directional selectivity')
            ax0.set_xlabel(r'$\theta$')
            ax0.set_ylabel('spikerate')

            ax1.set_title('Residuals')
            ax1.set_xlabel(r'$\theta$')
            ax1.set_ylabel('residuals')
            plt.annotate('SSE: '+str(sse[n]), xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')

    
    # Plot decay of reconstruction error
    if False:
        fig, ax = plt.subplots()
        ax.plot(numComps, err, label='MSE')
        ax.set_xlabel('Number of components')
        ax.set_ylabel(r'$||A - \hat{A}||^2$')
        ax.set_title('Decay of reconstruction error with number of principal components')
        now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        fig.savefig(plotDir+'subsetbMain_errPlot_'+now+'.pdf')

    plt.legend(loc='upper right')
    plt.show()











