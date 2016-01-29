#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: athul
# @Date:   2016-01-02 22:49:02
# @Last Modified by:   Athul
# @Last Modified time: 2016-01-17 14:07:59

from __future__ import division
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
plt.style.use('ggplot')

dataRoot = '../datasets/video/'
data = scipy.io.loadmat('../datasets/video/2013-28-06/1/AmpMov.mat')
data = data['AmpMov']

NumMovies       = data[0, 0]['NumMovies']
NumNeurons      = data[0, 0]['NumNeurons']
Sorted          = data[0, 0]['Sorted']
Blank           = data[0, 0]['Blank']

M_nat           = data[0, 0]['M_nat']
MT_nat          = data[0, 0]['MT_nat']
MTA_nat         = data[0, 0]['MTA_nat']
MTNA_nat        = data[0, 0]['MTNA_nat']

SC              = data[0, 0]['SC']
NC              = data[0, 0]['NC']

M_K0            = data[0, 0]['M_K0']
MT_K0           = data[0, 0]['MT_K0']
MTA_K0          = data[0, 0]['MTA_K0']
MTNA_K0         = data[0, 0]['MTNA_K0']

M_K1            = data[0, 0]['M_K1']
MT_K1           = data[0, 0]['MT_K1']
MTA_K1          = data[0, 0]['MTA_K1']
MTNA_K1         = data[0, 0]['MTNA_K1']

M_K1_5          = data[0, 0]['M_K1_5']
MT_K1_5         = data[0, 0]['MT_K1_5']
MTA_K1_5        = data[0, 0]['MTA_K1_5']
MTNA_K1_5       = data[0, 0]['MTNA_K1_5']

M_K2            = data[0, 0]['M_K2']
MT_K2           = data[0, 0]['MT_K2']
MTA_K2          = data[0, 0]['MTA_K2']
MTNA_K2         = data[0, 0]['MTNA_K2']
# ============================================
# Trials          = data[0, 0]['Trials']
# MP_nat          = data[0, 0]['MP_nat']
# T_nat           = data[0, 0]['T_nat']
# MP_K0           = data[0, 0]['MP_K0']
# TK0             = data[0, 0]['TK0']
# MP_K1           = data[0, 0]['MP_K1']
# TK1             = data[0, 0]['TK1']
# MP_K1_5         = data[0, 0]['MP_K1_5']
# TK1_5           = data[0, 0]['TK1_5']
# MP_K2           = data[0, 0]['MP_K2']
# TK2             = data[0, 0]['TK2']
# M_K3            = data[0, 0]['M_K3']
# MT_K3           = data[0, 0]['MT_K3']
# MTA_K3          = data[0, 0]['MTA_K3']
# MP_K3           = data[0, 0]['MP_K3']
# TK3             = data[0, 0]['TK3']
# CC              = data[0, 0]['CC']
# ED              = data[0, 0]['ED']
# MTNA_K3         = data[0, 0]['MTNA_K3']
# CC_SCNC         = data[0, 0]['CC_SCNC']
# =============================================