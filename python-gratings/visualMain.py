# -*- coding: utf-8 -*-
# @Author: Athul
# @Date:   2016-04-07 10:28:19
# @Last Modified by:   Athul
# @Last Modified time: 2016-04-07 11:45:49
from __future__ import division
import numpy as np
import scipy.io

import matplotlib.pyplot as plt
from matplotlib import animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Athul'), bitrate=1800)
plt.style.use('ggplot')

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

oriIdx = 1
matData = cellData[:60, oriIdx, :-1]
matData = np.reshape(matData, (6, 10, 120))

fig = plt.figure(figsize=(20, 12))
plt.axis('off')
plt.grid(False)
def f(a):
    global cellData
    return cellData[:, :, a]

im = plt.imshow(f(0), cmap=plt.get_cmap('gray'), interpolation='nearest', animated=True)

def updatefig(i):
    im.set_array(f(i))
    return im,
ani = animation.FuncAnimation(fig, updatefig, frames=np.arange(0, 121), interval=50, blit=True, repeat=False)
ani.save('ani.mp4', writer=writer)
plt.show()