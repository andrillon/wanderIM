#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 16:36:16 2018

@author: tand0009
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import kabuki

from patsy import dmatrix  # for generation of (regression) design matrices
from pandas import Series  # to manipulate data-frames generated by hddm


import hddm
print(hddm.__version__)

mydata = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_BehavPupil_ProbeResults.txt')
mydata.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stim','dprobe','pupil']
#mydata.columns = ['subj_idx','nblock','condition','ntrial','stimid','response','rt','stimulus','pupil']
#test.fillna(999, inplace=True)
mydata = mydata[np.isfinite(mydata['rt'])]


def z_link_func(x, data=mydata):
    stim = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
                               {'s': data.stim.ix[x.index]}))
    )
    return 1 / (1 + np.exp(-(x * stim)))

z_reg = {'model': 'z ~ 1 + C(pupil)', 'link_func': z_link_func}
v_reg = {'model': 'v ~ 1 + C(pupil)', 'link_func': lambda x: x}

reg_descr = [z_reg, v_reg]

m_reg = hddm.HDDMRegressor(mydata, reg_descr, include='z')

m_reg.sample(5000, burn=200)


mydata = mydata[np.logical_or(np.logical_and(mydata.dprobe > -17,mydata.stim ==0),np.logical_and(mydata.dprobe > -3,mydata.stim ==1))]
mydata2 = hddm.utils.flip_errors(mydata)
mydata2 = mydata2[mydata2.state != 4]
mydata2["state"] = mydata2["state"].astype('category')

m = hddm.models.HDDMRegressor(mydata2, 't ~ pupil')

m_within_subj = hddm.HDDMRegressor(mydata2, "t ~ C(stim, state('1'))")

m_reg = hddm.HDDMRegressor(mydata2,
                           "t ~ pupil:C(state)",
                           depends_on={'v': 'stim'})