# -*- coding: utf-8 -*-
# @Author: Athul Vijayan
# @Date:   2016-05-08 19:40:14
# @Last Modified by:   Athul Vijayan
# @Last Modified time: 2016-05-08 19:48:23
from __future__ import division
import numpy as np
import scipy.io
from lcs import rlcs as rlcs
import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

plotDir = '../plots/'

dataTargets = ['../datasets/driftingGratings/Mouse-A/', '../datasets/driftingGratings/Mouse-B/', '../datasets/driftingGratings/Mouse-C/', '../datasets/driftingGratings/Mouse-D/', '../datasets/driftingGratings/Mouse-E/']

mouse = 0
data = scipy.io.loadmat(dataTargets[mouse] + 'Data.mat')
data = data['Data']
smoothData = data[0, 0]['Spks']
stimuliSeq = data[0, 0]['StimSeq']

cellData = np.zeros((smoothData.shape[0], stimuliSeq.size, 121))
for i in xrange(smoothData.shape[0]):
    for j in xrange(stimuliSeq.size):
        cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
    # sort the array
    cellData[i] = cellData[i, np.argsort(cellData[i, :, -1])]
n1 = 10
X = cellData[n1, 10, :-1]

mouse = 3
data = scipy.io.loadmat(dataTargets[mouse] + 'Data-MouseD.mat')
data = data['Data']
smoothData = data[0, 0]['Spks']
stimuliSeq = data[0, 0]['StimSeq']

cellData = np.zeros((smoothData.shape[0], stimuliSeq.size, 121))
for i in xrange(smoothData.shape[0]):
    for j in xrange(stimuliSeq.size):
        cellData[i, j] = np.append(smoothData[i, 120*j: 120*(j+1)], np.radians(stimuliSeq[j]) )
    # sort the array
    cellData[i] = cellData[i, np.argsort(cellData[i, :, -1])]
n2 = 2
Y = cellData[n2, 10, :-1]

# ======================== RLCS start here ======================
tau_dist = 0.005
score, diag, cost = rlcs.rlcs(X, Y, tau_dist= tau_dist,  delta=0.5)

segment = rlcs.backtrack(X, Y, score, diag, cost)
lenSeg = segment.shape[0]

xSegs, ySegs = rlcs.getSoftSegments(segment, X, Y)
print len(xSegs)
# ========================= Plots here ===========================
# ================ Plot extracted subsequences ===================
fig, ax = plt.subplots()
lens = [i.shape[0] for i in xSegs]
idx = lens.index(max(lens))
xseg, yseg = xSegs[idx], ySegs[idx]
ax.plot(xrange(xseg.size), xseg)
ax.plot(xrange(yseg.size), yseg)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Matched signals across mice with dist_thres ' + str(tau_dist))
text = '''Neuron A = {0}\nNeuron B = {1}'''
ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
fig.savefig(plotDir+'rlcsMain_getSegs_mice'+now+'.eps')

# Plot the score matrix
fig, ax = rlcs.plotLCS(segment, X, Y)
ax.set_xlabel('template')
ax.set_ylabel('target')
ax.set_title('Match of signals across mice with dist_thres ' + str(tau_dist))
text = '''Neuron A = {0}\nNeuron B = {1}'''
ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
fig.savefig(plotDir+'rlcsMain_backtrack_mice_'+now+'.pdf')

plt.show()