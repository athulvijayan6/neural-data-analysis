#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-13 17:08:11
# @Last Modified by:   athul
# @Last Modified time: 2016-01-24 20:48:08
from __future__ import division
from scipy.stats import pearsonr
import numpy as np

def reliabilityCorr(spikeRates):
    corrCoef = 0
    count = 0
    for i in xrange(spikeRates.shape[1]-1):
        for j in xrange(i+1, spikeRates.shape[1]):
            corrCoef += pearsonr(spikeRates[:, i], spikeRates[:, j])[0]
            count += 1
    return corrCoef/count
    

if __name__ == '__main__':
    print 'reliability factor module'