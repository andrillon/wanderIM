#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 16:36:16 2018

@author: tand0009
"""

## Import required packages
import hddm
import numpy as np        
import pandas as pd
import matplotlib.pyplot as plt
from kabuki.analyze import gelman_rubin

# Load data from csv file into a NumPy structured array
data = mydata = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_BehavPupil_ProbeResults.txt')
data.columns = ['session','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stimulus','dprobe','pupil']

# Create histogram of RTs by subject, looks odd due to dummy coding of nogo response.
data = hddm.utils.flip_errors(data)
data = data[data.state != 4]
data["state"] = data["state"].astype('category')

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby(['task']):
    subj_data.rt.hist(bins=50, histtype='step', ax=ax)

#plt.savefig('hddm_RThists_bysubj.pdf')


# HDDMStimCoding Model Subclass enables you to specify your bounds as go/no-go, as opposed to correct/incorrect. This allows the estimation of both z bias and v bias. 
# Code below fits the full model where session (condition) influences all parameters, tends to have the lowest DIC. 
model = hddm.HDDMStimCoding(data, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'task', 'a': 'task', 't': 'task', 'z': 'task'}, p_outlier=.05)
model.find_starting_values()# Create model and start MCMC sampling
model.sample(5000, burn=500, thin=2, dbname='hddm_stimZ&v_gonogo.db', db='pickle')
model.save('hddm_stimZ&v_gonogo')
model.print_stats()


mydata = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_BehavPupil_ProbeResults.txt')
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