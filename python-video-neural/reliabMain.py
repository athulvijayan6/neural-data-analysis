#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-31 20:32:16
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-01 13:53:58
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
grid = {'centers': np.array([[ 1       ,  1 ], 
                             [ 2       ,  1 ],
                             [ 1.5     ,  1.5773 ]]),
              'x': np.array([ 5.]),
              'y': np.array([ 3.])}
d_matrix = np.random.rand(3)
ax = visualize.plot_map(grid, d_matrix)
# fig, ax = plt.subplots()
# centers = grid['centers']
# xpoints = centers[:, 0]
# ypoints = centers[:, 1]
# ax.scatter(xpoints, ypoints, color='red', marker='s')
# ax.axis([min(xpoints)-1., max(xpoints)+1., min(ypoints)-1., max(ypoints)+1.])

# xy_pixels = ax.transData.transform(centers)
# xpix, ypix = xy_pixels.T

# # discover radius and hexagon
# apothem = (xpix[1] - xpix[0]) / np.sqrt(3)
# area_inner_circle = np.pi * (apothem ** 2)
# collection_bg = RegularPolyCollection(
#     numsides=6,  # a hexagon
#     rotation=0,
#     sizes=(area_inner_circle,),
#     edgecolors = (0, 0, 0, 1),
#     array= d_matrix,
#     cmap = cm.Blues,
#     offsets = centers,
#     transOffset = ax.transData,
# )
# trans = transforms.Affine2D().scale(fig.dpi / 72.0)
# collection_bg.set_transform(trans)  # the points to pixels transform

# ax.add_collection(collection_bg, autolim=True)
# ax.autoscale_view()
# # plt.colorbar(collection_bg, cax=cax)



















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