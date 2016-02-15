#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-31 20:32:16
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-05 13:25:53
from __future__ import division
import datetime
import numpy as np
import scipy.io
from matplotlib import transforms
from reliability import reliability

from matplotlib import colors, cm
from matplotlib.collections import RegularPolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

from visualize import visualize
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']
plotDir = '../plots/'

NumNeurons      = data[0, 0]['NumNeurons']

MT_nat          = data[0, 0]['MT_nat']
MT_K0           = data[0, 0]['MT_K0']
MT_K1           = data[0, 0]['MT_K1']
MT_K1_5         = data[0, 0]['MT_K1_5']
MT_K3           = data[0, 0]['MT_K3']

vidIndex = 0
neuronId = 20
data = MT_nat
sample_rate = 20
ensembleSpikeRate = data[0, vidIndex]

# ========================== Reliablilty plot ================================
if True:
    reliables = []
    reliabilityIndex = []
    for neuronId in xrange(49):
        reliabililtyCorr = reliability.reliabilityCorr(ensembleSpikeRate[neuronId])
        reliabilityIndex.append(reliabililtyCorr)
        if reliabililtyCorr > 0.4:
            reliables.append(neuronId)
    reliabilityIndex = np.array(reliabilityIndex)
    fig, ax = plt.subplots()
    ax.scatter(xrange(reliabilityIndex.size), reliabilityIndex, color='red')
    ax.plot(xrange(reliabilityIndex.size), reliabilityIndex, color='black', alpha = 0.5)
    ax.set_xlabel('Neuron Id')
    ax.set_ylabel(r'$R_A$')
    ax.set_title('Reliability measure of neurons')

    now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    fig.savefig(plotDir+'reliabMain_raPlot_'+now+'.pdf')

# =============================================================================
# ============================ Reliability map ================================
plt.show() 