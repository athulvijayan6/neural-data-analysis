#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-02 22:49:02
# @Last Modified by:   athul
# @Last Modified time: 2016-01-25 15:57:15

from __future__ import division
import numpy as np
import scipy.io
from reliability import reliability
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']

NumMovies       = data[0, 0]['NumMovies']
NumNeurons      = data[0, 0]['NumNeurons']
Sorted          = data[0, 0]['Sorted']
Blank           = data[0, 0]['Blank']
Trials          = data[0, 0]['Trials']

# M_nat           = data[0, 0]['M_nat']
# M_K0            = data[0, 0]['M_K0']
# M_K1            = data[0, 0]['M_K1']
# M_K1_5          = data[0, 0]['M_K1_5']
# M_K2            = data[0, 0]['M_K2']
# M_K3            = data[0, 0]['M_K3']

MT_nat          = data[0, 0]['MT_nat']
MT_K0           = data[0, 0]['MT_K0']
MT_K1           = data[0, 0]['MT_K1']
MT_K1_5         = data[0, 0]['MT_K1_5']
MT_K3           = data[0, 0]['MT_K3']

# MTA_nat         = data[0, 0]['MTA_nat']
# MTA_K0          = data[0, 0]['MTA_K0']
# MTA_K1          = data[0, 0]['MTA_K1']
# MTA_K1_5        = data[0, 0]['MTA_K1_5']
# MTA_K2          = data[0, 0]['MTA_K2']
# MTA_K3          = data[0, 0]['MTA_K3']

# MTNA_nat        = data[0, 0]['MTNA_nat']
# MTNA_K0         = data[0, 0]['MTNA_K0']
# MTNA_K1         = data[0, 0]['MTNA_K1']
# MTNA_K1_5       = data[0, 0]['MTNA_K1_5']
# MTNA_K2         = data[0, 0]['MTNA_K2']
# MTNA_K3         = data[0, 0]['MTNA_K3']

# MP_nat          = data[0, 0]['MP_nat']
# MP_K0           = data[0, 0]['MP_K0']
# MP_K1           = data[0, 0]['MP_K1']
# MP_K1_5         = data[0, 0]['MP_K1_5']
# MP_K2           = data[0, 0]['MP_K2']
# MP_K3           = data[0, 0]['MP_K3']

# TK0             = data[0, 0]['TK0']   # Number of trials of each movie.
# TK1             = data[0, 0]['TK1']
# TK1_5           = data[0, 0]['TK1_5']
# TK2             = data[0, 0]['TK2']
# TK3             = data[0, 0]['TK3']

# CC              = data[0, 0]['CC']
# ED              = data[0, 0]['ED']
# CC_SCNC         = data[0, 0]['CC_SCNC']
# SC              = data[0, 0]['SC']
# NC              = data[0, 0]['NC']

# M_K0 --> 5 x 49 x 2000
# MT_K0 --> 5 x 49 x 200 x 10
# MTA_K0 --> 5 x 49 x 200
# MTNA_K0 --> 5 x 49 x 2000
# MP_K0 --> 5 x 9800 x 10
# ========================= Plot a single neuron response for same stimuli ==========
vidIndex = 0
neuronId = 20
data = MT_K0
sample_rate = 20
ensembleSpikeRate = data[0, vidIndex]

# ========================== Reliablilty plot ================================
reliables = []
reliabilityIndex = []
for neuronId in xrange(49):
    reliabililtyCorr = reliability.reliabilityCorr(ensembleSpikeRate[neuronId])
    reliabilityIndex.append(reliabililtyCorr)
    if reliabililtyCorr > 0.4:
        reliables.append(neuronId)
reliabilityIndex = np.array(reliabilityIndex)
fig, ax = plt.subplots()
heatmap = ax.pcolor(np.resize(reliabilityIndex, (1, reliabilityIndex.size)), cmap=plt.cm.Blues)
ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
print reliables

# ================================== response plots ======================
neuronId = 30
reliabililtyCorr = reliability.reliabilityCorr(ensembleSpikeRate[neuronId])
fig, ax = plt.subplots()
meanResponse = np.zeros(ensembleSpikeRate[0, :, 0].shape)
for trial in xrange(ensembleSpikeRate.shape[2]):
    spikeRate = ensembleSpikeRate[neuronId, :, trial]
    meanResponse += spikeRate
    x = [x/sample_rate for x in xrange(spikeRate.size)]
    x = np.array(x)
    ax.plot(x, spikeRate, label='trial '+str(trial+1), alpha=0.6)
meanResponse = meanResponse/(trial + 1)
ax.plot(x, meanResponse, label='Mean response', color='red')
ax.set_title('Responses to movie {0} for neuron {1} having stability correlation {2:.2f}'.format(vidIndex, neuronId, reliabililtyCorr))
ax.set_xlabel('time in seconds')
ax.set_ylabel('Spike rate')

# # =================================================================================
plt.legend(loc='upper right')
plt.show() 



