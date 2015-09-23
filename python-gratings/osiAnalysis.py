#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-04 16:42:24
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-23 12:40:54
from __future__ import division
import numpy as np
import scipy.io
import osi
from scipy.cluster.vq import kmeans, whiten, vq
import plotly.plotly as py
from plotly.graph_objs import *

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

        cirvar[i], dircirvar[i], OSI[i], DSI[i] = osi.computeAll(spikeRate[i])

     # % ----------------------- k-means clustering ----------
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
    
    showPlots = False
    if showPlots:
        # % ====================== Plots =======================
        neuronId = 63

        # % ----------------------- --------------------------------
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
        trace = []
        for i in range(3):
            trace.append(Scatter(x = oriData[i], y = dirData[i], mode='markers+text',
                name='Class '+str(i), marker=Marker(size=12,), text=names[i], textposition='top center'))
        plt = Data(trace)
        layout = Layout(title='Classification of simple and complex cells',
            xaxis=XAxis(title='Orientation index'),
            yaxis=YAxis(title='Direction index'))
        fig = Figure(data=plt, layout=layout)
        plot_url = py.plot(fig, filename='neuron-classification-kmeans')
        # ----------------------------------------------------------

        # % % 1. Plots each time-series data of length 120 with 
        # % % each trial of an experiment with seperate color
        theta = 1
        trace = []
        for j in xrange(10):
            trace.append(Scatter(x = list(xrange(1,121)), y = list(cellData[neuronId, 10*theta+j, :-2]), name='trial '+str(trial)))
        plt = Data(trace)
        layout = Layout(title='Plots of time-series data with each trial in diffenrent color',
            xaxis=XAxis(title='time'),
            yaxis=YAxis(title='smoothed amplitude'))
        fig = Figure(data=plt, layout=layout)
        plot_url = py.plot(fig, filename='neural-response-ts-frame')
        # # % --------------------------------------------------------

        # % --------------------------3-----------------------------
        # % 3. Plots response vs angle of stimulus pattern.
        trace = []
        for trial in xrange(10):
            trialData = spikeRate[neuronId, trial::10, :]
            trace.append(Scatter(x=list(trialData[:, 1]), y=list(trialData[:, 0]), name='trial '+str(trial)))
        plt = Data(trace)
        layout = Layout(title='Neuron response (spike rate) Vs direction of stimuli',
            xaxis=XAxis(title='angle in radians'),
            yaxis=YAxis(title='spike rate'))
        fig = Figure(data=plt, layout=layout)
        plot_url = py.plot(fig, filename='neural-response-vs-theta')
        # % --------------------------------------------------------


        # % --------------------------4-----------------------------
        # % 4. Polar plot of Directional selectivity
        # % Take a neuron, there will be 10 plots corresponding to each trial of the experiment
        # % plot each of the 10 plots in a figure to verify the similarity in response
        trace = []
        for trial in xrange(10):
            # % Normalize all trial response from 0 to 1
            normData = np.abs(spikeRate[neuronId, trial::10, 0])/np.max(np.abs(spikeRate[neuronId, trial::10, 0]));
            normData = np.append(normData, normData[0])
            theta = spikeRate[neuronId, trial::10, 1]
            theta = np.append(theta, theta[0])
            trace.append(Scatter(r=list(normData), t=list(theta), name='trial '+str(trial)))

        # draw dir vector
        plt = Data(trace)
        layout = Layout(title='Polar plot of Directional selectivity',)
        fig = Figure(data=plt, layout=layout)
        plot_url = py.plot(fig, filename='neural-response-dircirvar')

        # % --------------------------------------------------------

        # % --------------------------5-----------------------------
        # % 4. Polar plot of Orientation selectivity
        # % Take a neuron, there will be 10 plots corresponding to each trial of the experiment
        # % plot each of the 10 plots in a figure to verify the similarity in response

        # % There are observations for 16 angles, but there are only 8 orientations.
        # % We can treat them as two different observations. More data !

        trace = []
        for trial in xrange(10):
            # % Normalize all trial response from 0 to 1
            normData = np.abs(spikeRate[neuronId, trial::10, 0])/np.max(np.abs(spikeRate[neuronId, trial::10, 0]));
            normData = np.append(normData, normData[0])
            theta = 2*spikeRate[neuronId, trial::10, 1]
            theta = np.append(theta, theta[0])
            trace.append(Scatter(r=list(normData), t=list(theta), name='trial '+str(trial)))

        # draw ori vector
        plt = Data(trace)
        layout = Layout(title='Polar plot of Orientation selectivity',)
        fig = Figure(data=plt, layout=layout)
        plot_url = py.plot(fig, filename='neural-response-oricirvar')
        
        # % ********************** Plots ************************

