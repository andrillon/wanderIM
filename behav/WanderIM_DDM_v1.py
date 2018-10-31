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


import hddm
print(hddm.__version__)

data = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults.txt')
data.columns = ['subj_idx','nblock','task','ntrial','stimid','response','rt','stim']
data.fillna(999, inplace=True)
#data = data[np.isfinite(data['rt'])]
data2 = hddm.utils.flip_errors(data)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data2.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
    
    
# Instantiate model object passing it our data (no need to call flip_errors() before passing it).
# This will tailor an individual hierarchical DDM around your dataset.
m = hddm.HDDM(data2)
# find a good starting point which helps with the convergence.
m.find_starting_values()
# start drawing 7000 samples and discarding 5000 as burn-in
m.sample(2000, burn=20)   

m_stim = hddm.HDDM(data2, depends_on={'v': 'task'})
m_stim.find_starting_values()
m_stim.sample(5000, burn=3000)

kabuki.analyze.plot_posterior_nodes(m_stim.nodes_db.node[['v(1)', 'v(2)']])

probe = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults.txt')
probe.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stim','dprobe']
#probe.fillna(999, inplace=True)
#probe = probe[np.isfinite(probe['rt'])]
probe = probe[probe.state != 4]
probe = probe[np.logical_or(np.logical_and(probe.dprobe > -17,probe.stim ==0),np.logical_and(probe.dprobe > -3,probe.stim ==1))]
probe2 = hddm.utils.flip_errors(probe)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, state_data in probe2.groupby('state'):
    state_data.rt.hist(bins=50, histtype='step', ax=ax)
    
    
m_stim = hddm.HDDM(probe2, depends_on={'v': ['task','state'],'t': ['state']})
m_stim.find_starting_values()
m_stim.sample(5000, burn=3000)


kabuki.analyze.plot_posterior_nodes(m_stim.nodes_db.node[['v(1.2)', 'v(2.2)', 'v(3.2)']])
kabuki.analyze.plot_posterior_nodes(m_stim.nodes_db.node[['t(1)', 't(2)', 't(3)']])


probe = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults.txt')
probe.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stim','dprobe']
#probe.fillna(999, inplace=True)
probe = probe[np.isfinite(probe['rt'])]
probe = probe[probe.state != 4]
probe = probe[np.logical_or(np.logical_and(probe.dprobe > -17,probe.stim ==0),np.logical_and(probe.dprobe > -3,probe.stim ==1))]
probe2 = hddm.utils.flip_errors(probe)
  
    
m_stim2 = hddm.HDDM(probe2, depends_on={'v': ['task','state'],'t': ['state']})
m_stim2.find_starting_values()
m_stim2.sample(10000, burn=3000)


kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['v(1.1)', 'v(2.1)', 'v(3.1)']])
kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['v(1.2)', 'v(2.2)', 'v(3.2)']])
kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['t(1)', 't(2)', 't(3)']])

