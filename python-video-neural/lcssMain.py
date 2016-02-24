#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-23 16:00:11
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-24 16:00:57
from __future__ import division
import numpy as np
from lcs import lcss
import datetime
import scipy.io
from lcs import acf
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
X = np.mean(s1, axis=1)[:10]

Y = np.mean(s2, axis=1)[:10]


c, s, d, a, R = lcss.lcss(X, Y, thres_sim=0.4, thres_rc=0.3, rho=1)


# Plot the score matrix
if True:
    fig, ax = plt.subplots(figsize=(14, 12))
    cax = ax.imshow(d, aspect='auto', origin='lower', interpolation="none")
    ax.grid(True)
    cbar = fig.colorbar(cax)
    ax.set_xlabel('template')
    ax.set_ylabel('target')
    ax.set_title('Score matrix for pair of neurons')
    text = '''Neuron A = {0}\nNeuron B = {1}'''
    ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
    now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    ax.set_xticks(np.arange(12))
    ax.set_yticks(np.arange(12))
    # fig.savefig(plotDir+'lcssMain_scoremat_'+now+'.pdf')

plt.show()