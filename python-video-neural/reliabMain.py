#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-31 20:32:16
# @Last Modified by:   athul
# @Last Modified time: 2016-01-31 22:29:35
from __future__ import division
import numpy as np
import scipy.io
from reliability import reliability
from visualize import visualize
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']

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

# ========================== Reliablilty plot ================================
if False:
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

# =============================================================================
# ============================ Reliability map ================================
grid = {'centers': np.array([[ 1.5       ,  0.8660254 ],
   [ 2.5       ,  0.8660254 ],
   [ 3.5       ,  0.8660254 ],
   [ 4.5       ,  0.8660254 ],
   [ 5.5       ,  0.8660254 ],
   [ 6.5       ,  0.8660254 ],
   [ 1.        ,  1.73205081],
   [ 2.        ,  1.73205081],
   [ 3.        ,  1.73205081],
   [ 4.        ,  1.73205081],
   [ 5.        ,  1.73205081],
   [ 6.        ,  1.73205081],
   [ 1.5       ,  2.59807621],
   [ 2.5       ,  2.59807621],
   [ 3.5       ,  2.59807621],
   [ 4.5       ,  2.59807621],
   [ 5.5       ,  2.59807621],
   [ 6.5       ,  2.59807621],
   [ 1.        ,  3.46410162],
   [ 2.        ,  3.46410162],
   [ 3.        ,  3.46410162],
   [ 4.        ,  3.46410162],
   [ 5.        ,  3.46410162],
   [ 6.        ,  3.46410162],
   [ 1.5       ,  4.33012702],
   [ 2.5       ,  4.33012702],
   [ 3.5       ,  4.33012702],
   [ 4.5       ,  4.33012702],
   [ 5.5       ,  4.33012702],
   [ 6.5       ,  4.33012702],
   [ 1.        ,  5.19615242],
   [ 2.        ,  5.19615242],
   [ 3.        ,  5.19615242],
   [ 4.        ,  5.19615242],
   [ 5.        ,  5.19615242],
   [ 6.        ,  5.19615242]]),
'x': np.array([ 6.]),
'y': np.array([ 6.])}
d_matrix = np.random.rand(36)
ax = visualize.hexMap(grid, d_matrix, w=900)

reliables = []
reliabilityIndex = []
for neuronId in xrange(49):
    reliabililtyCorr = reliability.reliabilityCorr(ensembleSpikeRate[neuronId])
    reliabilityIndex.append(reliabililtyCorr)
    if reliabililtyCorr > 0.4:
        reliables.append(neuronId)
reliabilityIndex = np.array(reliabilityIndex)
# grid = {'centers':numNeurons, 'x': ,'y':}
print reliables

# plt.legend(loc='upper right')
plt.show() 