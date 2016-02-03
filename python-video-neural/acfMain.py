#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-02-03 11:36:02
# @Last Modified by:   Athul
# @Last Modified time: 2016-02-03 18:21:35

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

s1 = ensembleSpikeRate[0]
s2 = ensembleSpikeRate[1]

# average across trials
x = np.mean(s1, axis=1)
target = np.mean(s2, axis=1)


if False:
    width = 50
    nhat = width
    fig, ax = plt.subplots()
    template = acf.window(x, nhat, width)
    ccfunction, ax = acf.ccf(template, x, ax=ax)

width = 50
nhat = width
frame_shift = 10
while nhat <= target.size:
    template = acf.window(x, nhat, width)
    a = acf.ccf(template, target)
    try:
        acfGram = np.vstack((acfGram, a))
    except:
        acfGram = a
    nhat += frame_shift

plt.show() 