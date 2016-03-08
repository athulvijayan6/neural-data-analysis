#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-23 16:00:11
# @Last Modified by:   Athul
# @Last Modified time: 2016-03-08 12:21:41
from __future__ import division
import numpy as np
import scipy.io
from lcs import rlcs as rlcs
import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']
plotDir = '../plots/'

NumMovies       = data[0, 0]['NumMovies']
NumNeurons      = data[0, 0]['NumNeurons']

MT_nat          = data[0, 0]['MT_nat']
MT_K0           = data[0, 0]['MT_K0']
MT_K1           = data[0, 0]['MT_K1']
MT_K1_5         = data[0, 0]['MT_K1_5']
MT_K3           = data[0, 0]['MT_K3']

vidIndex = 0
neuronId = 20
data = MT_K0
sample_rate = 20
ensembleSpikeRate = data[0, vidIndex]
n1 , n2 = 10, 30
s1 = ensembleSpikeRate[n1]
s2 = ensembleSpikeRate[n2]

# average across trials
X = np.mean(s1, axis=1)
Y = np.mean(s2, axis=1)
# ======================== RLCS start here ======================
tau_dist = 0.005
score, diag, cost = rlcs.rlcs(X, Y, tau_dist= tau_dist,  delta=0.5)

segment = rlcs.backtrack(X, Y, score, diag, cost)
lenSeg = segment.shape[0]



# ========================= Plots here ===========================
# Plot the score matrix
if True:
    fig, ax = plt.subplots(figsize=(14, 12))
    cax = ax.imshow(score, aspect='auto', origin='lower', interpolation="none")
    ax.grid(True)
    cbar = fig.colorbar(cax)
    ax.set_xlabel('template')
    ax.set_ylabel('target')
    ax.set_title('Score matrix for pair of neurons  with dist_thres ' + str(tau_dist))
    text = '''Neuron A = {0}\nNeuron B = {1}'''
    ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
    now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    fig.savefig(plotDir+'rlcsMain_scoremat_'+now+'.pdf')

if True:
    match = np.zeros((X.size, Y.size))
    for m in segment:
        i, j = m[0], m[1]
        match[i-1, j-1] = m[2]
    fig, ax = plt.subplots(figsize=(14, 12))
    cax = ax.imshow(match, aspect='auto', origin='lower', interpolation="none")
    ax.grid(True)
    cbar = fig.colorbar(cax)
    ax.set_xlabel('template')
    ax.set_ylabel('target')
    ax.set_title('Match of signals after backtrack with dist_thres ' + str(tau_dist))
    text = '''Neuron A = {0}\nNeuron B = {1}'''
    ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
    now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    fig.savefig(plotDir+'rlcsMain_backtrack_'+now+'.pdf')


plt.show()