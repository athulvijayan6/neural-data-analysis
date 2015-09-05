#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2015-09-04 16:42:24
# @Last Modified by:   Athul
# @Last Modified time: 2015-09-05 19:10:09
from __future__ import division
import numpy as np
import scipy.io
import osi
from scipy.cluster.vq import kmeans, whiten, vq

dataTargets = ['dataset/Mouse-A/', 'dataset/Mouse-B/', 'dataset/Mouse-C/', 'dataset/Mouse-D/', 'dataset/Mouse-E/']
for mouse in xrange(5):
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
            cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
        # sort the array
        cellData[i] = cellData[i, np.argsort(cellData[i, :, -1])]
        # Calculate spike rate response of each neuron as a real number
        spikeRate[i] = osi.calculateSpikeRate(cellData[i])
        OSI_n = DSI_n = np.zeros((10))
        for trial in xrange(10):
            trialData = spikeRate[i, trial::10, :]
            R_pref_ori = np.max(trialData[:, 0])
            theta_pref_ori = np.argmax(trialData[:, 0])
            # compute OSI
            theta_orth = (theta_pref_ori + 4) % 16
            R_orth = trialData[theta_orth, 0]
            OSI_n[trial] = (R_pref_ori - R_orth)/(R_pref_ori + R_orth)
            # % Compute DSI
            theta_null = (theta_pref_ori + 8) % 16
            R_null = trialData[theta_null, 0]
            DSI_n[trial] = (R_pref_ori - R_null)/(R_pref_ori + R_null)
        OSI[i] = np.mean(OSI_n)
        DSI[i] = np.mean(DSI_n)
        # % compute orientation variance
        cirvar[i] = osi.cirVar(spikeRate[i])
        # % compute directional variance
        dircirvar[i] = osi.dirCirVar(spikeRate[i])

     # % ----------------------- k-means clustering ----------
    featureVec = np.vstack((OSI, DSI, cirvar, dircirvar))
    featureVec = np.transpose(featureVec)
    featureVec = whiten(featureVec)
    numClusters = 3
    centroids = kmeans(featureVec, numClusters)[0]
    idx = vq(featureVec, centroids)[0]

    # % Analyze the k means clustering results.
    globalMean = np.mean(featureVec, 0)
    total_ss = 0
    for i in xrange(featureVec.shape[0]):
        total_ss = total_ss + np.linalg.norm(globalMean - featureVec[i, :])
    between_ss = 0
    for i in xrange(numClusters):
        between_ss = between_ss + np.bincount(idx)[i]*np.linalg.norm(globalMean - centroids[i, :])
    clusterEfficiency = (between_ss/total_ss)*100
    print('The clustering is done with an efficiency of ' + str(clusterEfficiency))
    print('To see the definition of clustering efficiency, see the documentation')
    # % ----------------------- -----------------------------

    showPlots = False;
    if showPlots:
        # % ====================== Plots =======================
        colors = hsv(10);
        neuronId = 64;

        # % % 1. Plots each time-series data of length 120 with 
        # % % each trial of an experiment with seperate color
        # % figure;
        # % hold on;
        # % experimentNo = 2;
        # % for j=1:10
        # %     h = plot(1:120, cellData{1}(10*(experimentNo-1) + j, 1:120));
        # %     set(h,'Color', colors(j, :));
        # %     set(h,'LineWidth',2);
        # % end
        # % title('Plots of time-series data with each trial in diffenrent color');
        # % xlabel('time');
        # % ylabel('smoothed amplitude');
        # # % --------------------------------------------------------

        # % --------------------------3-----------------------------
        # % 3. Plots response vs angle of stimulus pattern.
        # figure;
        # for trial=4
        #     trialData = spikeRate{neuronId}(trial:10:end, :);
        #     h = plot(trialData(:, 2), trialData(:, 1));
        #     set(h,'Color', colors(trial-2, :));
        #     set(h,'Marker', 'o');
        #     set(h,'LineStyle', '--');
        #     set(h,'LineWidth', 2);
        #     hold on;
        # end
        # title('Neuron response (spike rate) Vs direction of stimuli');
        # ylabel('spike rate');
        # xlabel('angle in radians');
        # % --------------------------------------------------------


        # % --------------------------4-----------------------------
        # % 4. Polar plot of Directional selectivity
        # % Take a neuron, there will be 10 plots corresponding to each trial of the experiment
        # % plot each of the 10 plots in a figure to verify the similarity in response
        # figure;
        # for trial=1:10
        #     % Normalize all trial response from 0 to 1
        #     normData = abs(spikeRate{neuronId}(trial:10:end, 1))/max(abs(spikeRate{neuronId}(trial:10:end, 1)));
        #     normData(end+1) = normData(1);
        #     theta = spikeRate{neuronId}(trial:10:end, 2);
        #     theta(end+1) = theta(1);
        #     h = polar(theta, normData);
        #     set(h,'Color',colors(trial, :));
        #     set(h,'LineStyle','--');
        #     set(h,'Marker','o');
        #     set(h,'LineWidth',2);
        #     hold on;
        # end
        # h = compass(dircirvar(neuronId));
        # set(h,'LineWidth',3);
        # set(h,'Color', 'r');
        # title('Polar plot of Directional selectivity');
        # % --------------------------------------------------------

        # % --------------------------5-----------------------------
        # % 4. Polar plot of Orientation selectivity
        # % Take a neuron, there will be 10 plots corresponding to each trial of the experiment
        # % plot each of the 10 plots in a figure to verify the similarity in response

        # % There are observations for 16 angles, but there are only 8 orientations.
        # % We can treat them as two different observations. More data !

        # figure;
        # for trial=1:10
        #     % Normalize all trial response from 0 to 1
        #     normData = abs(spikeRate{neuronId}(trial:10:end, 1))/max(abs(spikeRate{neuronId}(trial:10:end, 1)));
        #     normData(end+1) = normData(1);
        #     theta = 2*spikeRate{neuronId}(trial:10:end, 2);
        #     theta(end+1) = theta(1);
        #     h = polar(theta, normData);
        #     set(h,'Color',colors(trial, :));
        #     set(h,'LineStyle','--');
        #     set(h,'Marker','o');
        #     set(h,'LineWidth',2);
        #     hold on;
        # end
        # h = compass(cirvar(neuronId));
        # set(h,'LineWidth',3);
        # set(h,'Color', 'r');
        # title('Polar plot of Orientation selectivity');
        # % --------------------------2-----------------------------

        # % ********************** Plots ************************

