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
#data = mydata = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_BehavPupil_ProbeResults.txt')
#data.columns = ['session','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','response','rt','stimulus','dprobe','pupil']
data = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults2.txt')
data.columns = ['subj_idx','nblock','task','ntrial','stimid','correctness','rt','stimulus','response','cond_v']
#data = data[np.logical_or(np.logical_and(data.dprobe > -17,data.stimulus ==0),np.logical_and(data.dprobe > -3,data.stimulus ==1))]
data.fillna(999, inplace=True)

# Create histogram of RTs by subject, looks odd due to dummy coding of nogo response.
data = hddm.utils.flip_errors(data)
#data = data[data.state != 4]
#data["state"] = data["state"].astype('category')

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby(['state']):
    subj_data.rt.hist(bins=50, histtype='step', ax=ax)

#plt.savefig('hddm_RThists_bysubj.pdf')


# HDDMStimCoding Model Subclass enables you to specify your bounds as go/no-go, as opposed to correct/incorrect. This allows the estimation of both z bias and v bias. 
# Code below fits the full model where session (condition) influences all parameters, tends to have the lowest DIC. 
model = hddm.HDDMStimCoding(data, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'cond_v', 'a': 'task', 't': 'task', 'z': 'task'}, p_outlier=.05)
model.find_starting_values()# Create model and start MCMC sampling
model.sample(500, burn=50, dbname='hddm_stim.db', db='pickle')
model.save('hddm_stim')
model.print_stats()

model.plot_posteriors(save=False)
model.plot_posterior_predictive()
model.plot_posteriors_conditions()

# Plot posterior probabilities for parameters and determine overlap
z_1, z_2 = model.nodes_db.node[['z(1)', 'z(2)']]
hddm.analyze.plot_posterior_nodes([z_1, z_2])
plt.xlabel('Response Bias z')
plt.ylabel('Posterior probability')
plt.title('Posterior of response bias between timepoints')
plt.savefig('hddm_demo_fig_06.pdf')

print("P(z(2) > z(1)) = ", (z_1.trace() > z_2.trace()).mean())

a_1, a_2 = model.nodes_db.node[['a(1)', 'a(2)']]
hddm.analyze.plot_posterior_nodes([a_1, a_2])
plt.xlabel('Boundary setting a')
plt.ylabel('Posterior probability')
plt.title('Posterior of boundary setting between timepoints')
plt.savefig('hddm_demo_fig_07.pdf')

print("P(a(2) > a(1)) = ", (a_1.trace() > a_2.trace()).mean())

t_1, t_2 = model.nodes_db.node[['t(1)', 't(2)']]
hddm.analyze.plot_posterior_nodes([t_1, t_2])
plt.xlabel('NDT')
plt.ylabel('Posterior probability')
plt.title('Posterior of NDT between timepoints')
plt.savefig('hddm_demo_fig_08.pdf')

print("P(t(2) < t(1)) = ", (t_1.trace() < t_2.trace()).mean())

v_1, v_2 = model.nodes_db.node[['v(goF)', 'v(goD)']]
hddm.analyze.plot_posterior_nodes([v_1, v_2])
plt.xlabel('Drift')
plt.ylabel('Posterior probability')
plt.title('Posterior of Drift between timepoints')
plt.savefig('hddm_demo_fig_09.pdf')

print("P(v(goD) > v(goF)) = ", (v_1.trace() > v_2.trace()).mean())