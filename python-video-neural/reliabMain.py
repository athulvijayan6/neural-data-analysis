#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-31 20:32:16
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-04 15:05:30
from __future__ import division
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
    heatmap = ax.pcolor(np.resize(reliabilityIndex, (1, reliabilityIndex.size)), cmap=plt.cm.Blues)
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    print reliables

# =============================================================================
# ============================ Reliability map ================================
# w = 10
# hits = np.array([1, 24, 14, 16,  6, 11,  8, 23, 15, 16, 15,  9, 20,  1,  3, 29,  4,
#               32, 22,  7, 26, 26, 35, 23,  7,  6, 11,  9, 18, 17, 22, 19, 34,  1,
#               36,  3, 31, 10, 22, 11, 21, 18, 29,  3,  6, 32, 15, 30, 27])
# n_centers = np.array([[ 1.5       ,  0.8660254 ],
#                  [ 2.5       ,  0.8660254 ],
#                  [ 3.5       ,  0.8660254 ],
#                  [ 4.5       ,  0.8660254 ],
#                  [ 5.5       ,  0.8660254 ],
#                  [ 6.5       ,  0.8660254 ],
#                  [ 1.        ,  1.73205081],
#                  [ 2.        ,  1.73205081],
#                  [ 3.        ,  1.73205081],
#                  [ 4.        ,  1.73205081],
#                  [ 5.        ,  1.73205081],
#                  [ 6.        ,  1.73205081],
#                  [ 1.5       ,  2.59807621],
#                  [ 2.5       ,  2.59807621],
#                  [ 3.5       ,  2.59807621],
#                  [ 4.5       ,  2.59807621],
#                  [ 5.5       ,  2.59807621],
#                  [ 6.5       ,  2.59807621],
#                  [ 1.        ,  3.46410162],
#                  [ 2.        ,  3.46410162],
#                  [ 3.        ,  3.46410162],
#                  [ 4.        ,  3.46410162],
#                  [ 5.        ,  3.46410162],
#                  [ 6.        ,  3.46410162],
#                  [ 1.5       ,  4.33012702],
#                  [ 2.5       ,  4.33012702],
#                  [ 3.5       ,  4.33012702],
#                  [ 4.5       ,  4.33012702],
#                  [ 5.5       ,  4.33012702],
#                  [ 6.5       ,  4.33012702],
#                  [ 1.        ,  5.19615242],
#                  [ 2.        ,  5.19615242],
#                  [ 3.        ,  5.19615242],
#                  [ 4.        ,  5.19615242],
#                  [ 5.        ,  5.19615242],
#                  [ 6.        ,  5.19615242]])

# fig, ax = plt.subplots()

# for p in n_centers:
#   ax.add_patch(
#     patches.RegularPolygon(
#         (p[0], p[1]),     # (x,y)
#         6,              # number of vertices
#         0.2,            # radius
#     )
#   )
















# reliables = []
# reliabilityIndex = []
# for neuronId in xrange(49):
#     reliabililtyCorr = reliability.reliabilityCorr(ensembleSpikeRate[neuronId])
#     reliabilityIndex.append(reliabililtyCorr)
#     if reliabililtyCorr > 0.4:
#         reliables.append(neuronId)
# reliabilityIndex = np.array(reliabilityIndex)
# # grid = {'centers':numNeurons, 'x': ,'y':}
# print reliables

# plt.legend(loc='upper right')
plt.show() 