# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-05-05 14:51:06
# @Last Modified by:   Athul
# @Last Modified time: 2016-05-05 14:55:00
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-23 16:00:11
# @Last Modified by:   Athul
# @Last Modified time: 2016-05-05 14:50:04
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
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']


NumMovies       = data[0, 0]['NumMovies']
NumNeurons      = data[0, 0]['NumNeurons']

MT_nat          = data[0, 0]['MT_nat']

vidIndex = 0
neuronId = 20
data = MT_nat
sample_rate = 20
ensembleSpikeRate = data[0, vidIndex]
n1 , n2 = 10, 30
s1 = ensembleSpikeRate[n1]
s2 = ensembleSpikeRate[n2]

# average across trials
X = np.mean(s1, axis=1)
Y = np.mean(s2, axis=1)

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

ax.fill_between(xrange(xseg.size), yseg, xseg, where=(yseg - xseg)**2 < 2, facecolor='green', interpolate=True)
plt.show()