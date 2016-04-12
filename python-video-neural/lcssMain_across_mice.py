#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-23 16:00:11
# @Last Modified by:   Athul Vijayan
# @Last Modified time: 2016-04-12 08:13:51
from __future__ import division
import numpy as np
import scipy.io
from lcs import rlcs as rlcs
import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

plotDir = '../plots/'

# ============================ Loading neuronal data here ===============
data = scipy.io.loadmat('../datasets/video/2014-03-04/1/AmpMov.mat')
data = data['AmpMov']

NumNeurons      = data[0, 0]['NumNeurons']

MT_nat          = data[0, 0]['MT_nat']

vidIndex = 0
data = MT_nat
ensembleSpikeRate = data[0, vidIndex]
n1 = 30
s = ensembleSpikeRate[n1]
# average across trials
X = np.mean(s, axis=1)

data = scipy.io.loadmat('../datasets/video/2013-23-06/2/AmpMov.mat')
data = data['AmpMov']
NumNeurons      = data[0, 0]['NumNeurons']
MT_nat          = data[0, 0]['MT_nat']
vidIndex = 0
data = MT_nat
ensembleSpikeRate = data[0, vidIndex]
n2 = 30
s = ensembleSpikeRate[n2]
# average across trials
Y = np.mean(s, axis=1)



# At the end of loading your data, Have query as X and reference as Y
# Both X and Y are numpy arrays
# In 1D case, both will be vectors
# in multidimensional, rows of X and Y denote each sample.
# and columns denote feature dimension
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
fig.savefig(plotDir+'rlcsMain_getSegs_'+now+'.eps')

# Plot the score matrix
fig, ax = rlcs.plotLCS(segment, X, Y)
ax.set_xlabel('template')
ax.set_ylabel('target')
ax.set_title('Match of signals across mice with dist_thres ' + str(tau_dist))
text = '''Neuron A = {0}\nNeuron B = {1}'''
ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
fig.savefig(plotDir+'rlcsMain_backtrack_'+now+'.pdf')

plt.show()