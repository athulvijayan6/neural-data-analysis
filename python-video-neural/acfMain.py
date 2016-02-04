#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-03 11:36:02
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-04 15:33:09

from __future__ import division
import datetime
import numpy as np
import scipy.io
from lcs import acf
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']
plotDir = '../plots/'

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
n1 , n2 = 0, 0
s1 = ensembleSpikeRate[n1]
s2 = ensembleSpikeRate[n2]

# average across trials
x = np.mean(s1, axis=1)
target = np.mean(s2, axis=1)

width = 140
nhat = width
if True:
    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax1 = plt.subplot(gs[0])
    template = acf.window(x, nhat, width)
    ccfunction, ax1 = acf.ccf(template, target, ax=ax1)
    text = '''Template neuron = {0}\nTarget neuron = {1}\nframe width = {2}\nframe ending = {3}'''
    ax1.set_title('Cross correlation of template and target')
    ax1.annotate(text.format(n1, n2, width, nhat), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
    ax2 = plt.subplot(gs[1])
    t = np.arange(template.size)
    ax2.plot(t, template)
    ax2.set_xlabel('samples')
    ax2.set_ylabel('Spike rate')
    ax2.set_title('Template signal')
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.savefig(plotDir+'acfMain_acfPlot_'+now+'.pdf')

if True:
    frame_shift = 5
    while nhat <= target.size:
        template = acf.window(x, nhat, width)
        a = acf.ccf(template, target)
        try:
            acfGram = np.vstack((acfGram, a))
        except:
            acfGram = a
        nhat += frame_shift
    acfGram = np.transpose(acfGram)
    fig, ax = plt.subplots(figsize=(14, 10))
    cax = ax.imshow(acfGram, aspect='auto', origin='lower')
    ax.grid(True)
    cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
    ax.set_xlabel('Frame ending (samples)')
    ax.set_ylabel('Lag')
    text = '''Neuron A = {0}\nNeuron B = {1}\nframe width = {2}\nframe shift = {3}'''
    ax.set_title('Correlationgram of pair of neurons')
    ax.annotate(text.format(n1, n2, width, frame_shift), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.savefig(plotDir+'acfMain_corrGram_'+now+'.eps')

plt.show() 