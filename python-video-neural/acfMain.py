#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-03 11:36:02
# @Last Modified by:   athul
# @Last Modified time: 2016-02-04 09:57:47

from __future__ import division
import numpy as np
import scipy.io
from lcs import acf
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']

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
n1 , n2 = 0, 1
s1 = ensembleSpikeRate[n1]
s2 = ensembleSpikeRate[n2]

# average across trials
x = np.mean(s1, axis=1)
target = np.mean(s2, axis=1)

width = 20
nhat = width
if True:
    fig, ax = plt.subplots()
    template = acf.window(x, nhat, width)
    ccfunction, ax = acf.ccf(template, x, ax=ax)

frame_shift = 10
while nhat <= target.size:
    template = acf.window(x, nhat, width)
    a = acf.ccf(template, x)
    try:
        acfGram = np.vstack((acfGram, a))
    except:
        acfGram = a
    nhat += frame_shift
acfGram = np.transpose(acfGram)
fig, ax = plt.subplots()
cax = ax.imshow(acfGram, aspect='auto', origin='lower', interpolation='none')
ax.grid(True)
cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
ax.set_xlabel('time')
ax.set_ylabel('Correlation function')
ax.set_title('Correlationgram of neuron '+str(n1)+' and neuron '+str(n2))
plt.show() 