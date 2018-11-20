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

cd /Users/tand0009/WorkGit/projects/inprogress/wanderIM/behav

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
probe.fillna(999, inplace=True)
#probe = probe[np.isfinite(probe['rt'])]
probe = probe[probe.state != 4]
probe = probe[np.logical_or(np.logical_and(probe.dprobe > -17,probe.stim ==0),np.logical_and(probe.dprobe > -3,probe.stim ==1))]
probe2 = hddm.utils.flip_errors(probe)
probe2['state'] = probe2['state'].astype('category')
probe2['task'] = probe2['task'].astype('category')
probe2['stim'] = probe2['stim'].astype('category')
                             
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, state_data in probe2.groupby('state'):
    state_data.rt.hist(bins=50, histtype='step', ax=ax)
    
    
    
m_stim = hddm.HDDM(probe2, depends_on={'v': ['task','state'],'t': ['state']})
m_stim.find_starting_values()
m_stim.sample(5000, burn=2000, dbname='traces.db', db='pickle')


kabuki.analyze.plot_posterior_nodes(m_stim.nodes_db.node[['v(1.1)', 'v(2.1)', 'v(3.1)']])
kabuki.analyze.plot_posterior_nodes(m_stim.nodes_db.node[['v(1.2)', 'v(2.2)', 'v(3.2)']])
kabuki.analyze.plot_posterior_nodes(m_stim.nodes_db.node[['t(1)', 't(2)', 't(3)']])

#probe = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults.txt')
#probe.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stim','dprobe']
#probe.fillna(999, inplace=True)
#probe = probe[np.isfinite(probe['rt'])]
#probe = probe[probe.state != 4]
#probe = probe[np.logical_or(np.logical_and(probe.dprobe > -17,probe.stim ==0),np.logical_and(probe.dprobe > -3,probe.stim ==1))]
#probe2 = hddm.utils.flip_errors(probe)
  
    
m_stim2 = hddm.HDDM(probe2, depends_on={'v': ['task','state'],'t': ['task','state']})
m_stim2.find_starting_values()
m_stim2.sample(5000, burn=1000, dbname='traces.db', db='pickle')
m_stim2.save('wIM_hddm_model_vt_taskstate2')


kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['v(1.1)', 'v(2.1)', 'v(3.1)']])
plt.savefig('hddm_model_vt_taskstate_v_F.pdf')
kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['v(1.2)', 'v(2.2)', 'v(3.2)']])
plt.savefig('hddm_model_vt_taskstate_v_D.pdf')  
kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['t(1.1)', 't(2.1)', 't(3.1)']])
plt.savefig('hddm_model_vt_taskstate_t_F.pdf')
kabuki.analyze.plot_posterior_nodes(m_stim2.nodes_db.node[['t(1.2)', 't(2.2)', 't(3.2)']])
plt.savefig('hddm_model_vt_taskstate_t_D.pdf')

v_ON, v_MW, v_MB= m_stim2.nodes_db.node[['v(1.1)', 'v(2.1)', 'v(3.1)']]
print('v-Face: (1) ON vs MW, (2) ON vs MB, (3) MW vs MB')
print((v_ON.trace() > v_MW.trace()).mean())
print((v_ON.trace() > v_MB.trace()).mean())
print((v_MW.trace() > v_MB.trace()).mean())
v_ON, v_MW, v_MB= m_stim2.nodes_db.node[['v(1.2)', 'v(2.2)', 'v(3.2)']]
print('v-Digit: (1) ON vs MW, (2) ON vs MB, (3) MW vs MB')
print((v_ON.trace() > v_MW.trace()).mean())
print((v_ON.trace() > v_MB.trace()).mean())
print((v_MW.trace() > v_MB.trace()).mean())

t_ON, t_MW, t_MB= m_stim2.nodes_db.node[['v(1.1)', 'v(2.1)', 'v(3.1)']]
print('t-Face: (1) ON vs MW, (2) ON vs MB, (3) MW vs MB')
print((t_ON.trace() > t_MW.trace()).mean())
print((t_ON.trace() > t_MB.trace()).mean())
print((t_MW.trace() > t_MB.trace()).mean())
t_ON, t_MW, t_MB= m_stim2.nodes_db.node[['v(1.2)', 'v(2.2)', 'v(3.2)']]
print('t-Digit: (1) ON vs MW, (2) ON vs MB, (3) MW vs MB')
print((t_ON.trace() > t_MW.trace()).mean())
print((t_ON.trace() > t_MB.trace()).mean())
print((t_MW.trace() > t_MB.trace()).mean())

from patsy import dmatrix
dmatrix("C(state, Treatment(1))", probe2.head(10))

probe = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults.txt')
probe.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stim','dprobe']
probe = probe[np.isfinite(probe['rt'])]
probe = probe[probe.state != 4]
probe = probe[np.logical_or(np.logical_and(probe.dprobe > -17,probe.stim ==0),np.logical_and(probe.dprobe > -3,probe.stim ==1))]
probe3 = hddm.utils.flip_errors(probe)
probe3['state'] = probe3['state'].astype('category')
probe3['task'] = probe3['task'].astype('category')
probe3['stim'] = probe3['stim'].astype('category')

m_within_subj = hddm.HDDMRegressor(probe3, "v ~ C(state, Treatment(1))")
m_within_subj.sample(5000, burn=1000)

v_ON, v_MW, v_MB = m_within_subj.nodes_db.ix[["v_Intercept",
                                              "v_C(state, Treatment(1))[T.2]",
                                              "v_C(state, Treatment(1))[T.3]"], 'node']
hddm.analyze.plot_posterior_nodes([v_ON, v_MW, v_MB])
print((v_ON.trace() > v_MW.trace()).mean())
print((v_ON.trace() > v_MB.trace()).mean())
print((v_MW.trace() > 0).mean())
print((v_MB.trace() > 0).mean())

m_within_subj2 = hddm.HDDMRegressor(probe3, "v ~ C(state, Treatment(1))", depends_on={'v': ['task','state'],'t': ['task','state']})
m_within_subj2.sample(5000, burn=1000)
