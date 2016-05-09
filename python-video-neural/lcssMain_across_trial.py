#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-23 16:00:11
# @Last Modified by:   Athul Vijayan
# @Last Modified time: 2016-05-08 19:18:25
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
dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2014-08-04/1/AmpMov.mat')
data = data['AmpMov']


NumMovies       = data[0, 0]['NumMovies']
NumNeurons      = data[0, 0]['NumNeurons']

MT_nat          = data[0, 0]['MT_nat']

vidIndex = 0
neuronId = 20
data = MT_nat
sample_rate = 20
ensembleSpikeRate = data[0, vidIndex]
# n = 10
n = 42
s = ensembleSpikeRate[n]
# average across trials
X = s[:, 0]
Y = s[:, 5]

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
ax.set_title('Matched signals for rlcs with dist_thres ' + str(tau_dist))
now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
fig.savefig(plotDir+'rlcsMain_motif_trial_0_1n_'+str(n)+now+'.eps')

# Plot the score matrix
fig, ax = rlcs.plotLCS(segment, X, Y)
ax.set_xlabel('template')
ax.set_ylabel('target')
ax.set_title('Score matrix from RLCS with dist_thres ' + str(tau_dist))
now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
fig.savefig(plotDir+'rlcsMain_score_trial_0_1n_'+str(n)+now+'.pdf')

plt.show()