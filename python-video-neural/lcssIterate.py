# -*- coding: utf-8 -*-
# @Author: Athul Vijayan
# @Date:   2016-04-05 23:39:31
# @Last Modified by:   Athul Vijayan
# @Last Modified time: 2016-04-12 08:28:15
from __future__ import division
import numpy as np
import scipy.io
import itertools
import re, os
import pickle
from lcs import rlcs as rlcs
import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')
plt.ioff()

dataRoot = '../datasets/'

import glob

datapaths = glob.glob('../datasets/**/**/**/AmpMov.mat')

tau_dist = 0.005

allSegs = {}
for path in datapaths:

    plotDir = '../plots/'+re.search("[0-9\-]+", path).group()+'/'
    print plotDir
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)

    data    = scipy.io.loadmat(path)
    data    = data['AmpMov']
    NumNeurons = data[0, 0]['NumNeurons']
    MT_nat  = data[0, 0]['MT_nat']
    
    vidIndex = 0
    data    = MT_nat
    sample_rate = 20
    ensembleSpikeRate = data[0, vidIndex]

    miceSegs = {}
    for n1, n2 in [j for j in itertools.combinations(xrange(NumNeurons), 2)]:
        s1 = ensembleSpikeRate[n1]
        s2 = ensembleSpikeRate[n2]

        # average across trials
        X = np.mean(s1, axis=1)
        Y = np.mean(s2, axis=1)
        # ======================== RLCS start here ======================
        score, diag, cost = rlcs.rlcs(X, Y, tau_dist= tau_dist,  delta=0.5)
        segment = rlcs.backtrack(X, Y, score, diag, cost)
        xSegs, ySegs = rlcs.getSoftSegments(segment, X, Y)
        # ========================= Plots here ===========================
        # ================ Plot extracted subsequences ===================

        fig, ax = plt.subplots()
        lens = [i.shape[0] for i in xSegs]
        idx = lens.index(max(lens))
        xseg, yseg = xSegs[idx], ySegs[idx]
        miceSegs[(n1, n2)] = xSegs
        
        ax.plot(xrange(xseg.size), xseg)
        ax.plot(xrange(yseg.size), yseg)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Matched signals for rlcs with dist_thres ' + str(tau_dist))
        text = '''Neuron A = {0}\nNeuron B = {1}'''
        ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
        now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        fig.savefig(plotDir+'rlcsMain_getSegs_'+now+'.eps')
        plt.close()

        # Plot the score matrix
        fig, ax = rlcs.plotLCS(segment, X, Y)
        ax.set_xlabel('template')
        ax.set_ylabel('target')
        ax.set_title('Match of signals after backtrack with dist_thres ' + str(tau_dist))
        text = '''Neuron A = {0}\nNeuron B = {1}'''
        ax.annotate(text.format(n1, n2), xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)
        now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        fig.savefig(plotDir+'rlcsMain_backtrack_'+now+'.pdf')
        plt.close()

    allSegs[path] = miceSegs
with open('outfile.pkl', 'wb') as f:
    pickle.dump(allSegs, f)

plt.show()
