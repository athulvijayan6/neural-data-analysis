# -*- coding: utf-8 -*-
# @Author: Athul Vijayan
# @Date:   2016-05-04 06:48:15
# @Last Modified by:   Athul Vijayan
# @Last Modified time: 2016-05-04 07:57:33
from __future__ import division
import datetime
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

data = scipy.io.loadmat('../../datasets/video/2013-28-06/2/AmpMov.mat')
data = data['AmpMov']
plotDir = '../plots/'

MT_nat = data[0, 0]['MT_nat']
vidIndex = 0
neuronId = 20
data = MT_nat
sample_rate = 20
ensembleSpikeRate = data[0, vidIndex]
n1 , n2 = 25, 15
s1 = ensembleSpikeRate[n1]
s2 = ensembleSpikeRate[n2]

# average across trials
x = np.mean(s1, axis=1)
target = np.mean(s2, axis=1)

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True, figsize=(10, 4))
axarr[0].plot(xrange(x.size), target)
axarr[0].set_title('Target sequence')
axarr[0].set_ylabel('Spike rate')
axarr[1].plot(xrange(x.size), x, alpha=0.6)
nhat = 80
width = 40
axarr[1].set_title('Extracting search frame from template')
axarr[1].set_ylabel('Spike rate')
axarr[1].plot(xrange(nhat-width,nhat), x[nhat-width:nhat], linewidth=1)
axarr[1].annotate('Ending of selected frame', xy=(nhat, x[nhat]), xycoords='data',
                xytext=(0.8, 0.95), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.05),
                horizontalalignment='right', verticalalignment='top', fontsize=15
                )
axarr[0].set_axis_bgcolor('white')
axarr[1].set_axis_bgcolor('white')
axarr[0].grid(True)
axarr[1].grid(True)
plt.show()
