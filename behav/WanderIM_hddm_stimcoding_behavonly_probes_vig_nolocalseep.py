#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 16:36:16 2018

@author: tand0009
"""

## Import required packages
import hddm
import numpy as np        
import scipy as sp        
import pandas as pd
import matplotlib.pyplot as plt
from kabuki.analyze import gelman_rubin


# Load data from csv file into a NumPy structured array
data = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults2.txt')
data.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','corr','rt','stim','dprobe','response','stimulus','vigC','cond_v']
#data = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults2.txt')
#data.columns = ['subj_idx','nblock','task','ntrial','stimid','correctness','rt','stimulus','response','cond_v']
data = data[np.logical_or(np.logical_and(data.dprobe > -19,data.stimulus ==1),np.logical_and(data.dprobe > -3,data.stimulus ==0))]
data.fillna(999, inplace=True)
data = data[data.rt > 0.2]

# Create histogram of RTs by subject, looks odd due to dummy coding of nogo response.
data = hddm.utils.flip_errors(data)
data = data[data.state != 4]
data["state"] = data["state"].astype('category')



#plt.savefig('hddm_RThists_bysubj.pdf')

dataD = data[data.task == 2]
dataF = data[data.task == 1]

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dataD.groupby(['vigC']):
    subj_data.rt.hist(bins=25, histtype='step', ax=ax)
    
# HDDMStimCoding Model Subclass enables you to specify your bounds as go/no-go, as opposed to correct/incorrect. This allows the estimation of both z bias and v bias. 
# Code below fits the full model where session (condition) influences all parameters, tends to have the lowest DIC. 
modelD = hddm.HDDMStimCoding(dataD, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'cond_v', 'a': 'vigC', 't': 'vigC', 'z': 'vigC'}, p_outlier=.05)
modelD.find_starting_values()# Create model and start MCMC sampling
modelD.sample(10000, burn=1000, dbname='hddm_stim.db', db='pickle')
#modelD.save('hddm_stim')
modelD.print_stats()

modelD.plot_posteriors(save=False)
modelD.plot_posterior_predictive()
modelD.plot_posteriors_conditions()

# Plot posterior probabilities for parameters and determine overlap
z_1, z_2 = modelD.nodes_db.node[['z(0)', 'z(1)']]
hddm.analyze.plot_posterior_nodes([z_1, z_2])
plt.xlabel('Response Bias z')
plt.ylabel('Posterior probability')
plt.title('Posterior of response bias between timepoints')

print("P(z(2) > z(1)) = ", (z_1.trace() > z_2.trace()).mean())

a_1, a_2 = modelD.nodes_db.node[['a(0)', 'a(1)']]
hddm.analyze.plot_posterior_nodes([a_1, a_2])
plt.xlabel('Boundary setting a')
plt.ylabel('Posterior probability')
plt.title('Posterior of boundary setting between timepoints')

print("P(a(2) > a(1)) = ", (a_1.trace() > a_2.trace()).mean())

t_1, t_2 = modelD.nodes_db.node[['t(0)', 't(1)']]
hddm.analyze.plot_posterior_nodes([t_1, t_2])
plt.xlabel('NDT')
plt.ylabel('Posterior probability')
plt.title('Posterior of NDT between timepoints')

print("P(t(2) < t(1)) = ", (t_1.trace() < t_2.trace()).mean())


v_g1, v_g2, v_ng1, v_ng2 = modelD.nodes_db.node[['v(go_alert)', 'v(go_drowsy)','v(nogo_alert)', 'v(nogo_drowsy)']]
hddm.analyze.plot_posterior_nodes([v_g1, v_g2, v_ng1, v_ng2])
plt.xlabel('Drift')
plt.ylabel('Posterior probability')
plt.title('Posterior of Drift between timepoints')

print("P(v(go_alert) > v(go_drowsy)) = ", (v_g1.trace() > v_g2.trace()).mean())
print("P(v(nogo_alert) > v(nogo_drowsy)) = ", (v_ng1.trace() > v_ng2.trace()).mean())

plt.plot((v_g1.trace() + v_ng1.trace()))
plt.plot((v_g2.trace() + v_ng2.trace()))

fig = plt.figure()
dc1=(v_g1.trace() + v_ng1.trace())
dc2=(v_g2.trace() + v_ng2.trace())
plt.hist(dc1,color='blue')
plt.hist(dc2,color='red')

mean_v_g1=(v_g2.trace()).mean()
mean_v_g0=(v_g1.trace()).mean()
mean_v_ng1=(v_ng2.trace()).mean()
mean_v_ng0=(v_ng1.trace()).mean()
mean_z_1=(z_2.trace()).mean()
mean_z_0=(z_1.trace()).mean()
mean_a_1=(a_2.trace()).mean()
mean_a_0=(a_1.trace()).mean()
mean_t_1=(t_2.trace()).mean()
mean_t_0=(t_1.trace()).mean()

std_v_g1=(v_g2.trace()).std()
std_v_g0=(v_g1.trace()).std()
std_v_ng1=(v_ng2.trace()).std()
std_v_ng0=(v_ng1.trace()).std()
std_z_1=(z_2.trace()).std()
std_z_0=(z_1.trace()).std()
std_a_1=(a_2.trace()).std()
std_a_0=(a_1.trace()).std()
std_t_1=(t_2.trace()).std()
std_t_0=(t_1.trace()).std()

x_pos=np.array([1,1,2,2,3,3,4,4,5,5])+np.array([-0.2,0.2,-0.2,0.2,-0.2,0.2,-0.2,0.2,-0.2,0.2])
xtick_pos=np.array([1,2,3,4,5])
xtick_labels = ['v_{go}', 'v_{nogo}', 'z', 'a', 't']
means = [mean_v_g0, mean_v_g1, mean_v_ng0, mean_v_ng1, mean_z_0, mean_z_1, mean_a_0, mean_a_1, mean_t_0, mean_t_1]
errors = [std_v_g0, std_v_g1, std_v_ng0, std_v_ng1, std_z_0, std_z_1, std_a_0, std_a_1, std_t_0, std_t_1]

fig, ax = plt.subplots()
ax.bar(x_pos, means, yerr=errors, align='center', alpha=0.5, ecolor='black', capsize=10, width=0.38, color=['blue','red','blue','red','blue','red','blue','red','blue','red'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos)
ax.set_xticklabels(xtick_labels)
ax.axhline(y=0,color='black',linewidth=0.5)
ax.set_title('Subj Rating informed HDDM')
ax.yaxis.grid(True)

fig, ax = plt.subplots()
ax.bar(x_pos[0:4], means[0:4], yerr=errors[0:4], align='center', alpha=0.5, ecolor='black', capsize=10, width=0.38, color=['blue','red','blue','red'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos[0:2])
ax.set_xticklabels(xtick_labels[0:2])
ax.set_title('Subj Rating informed HDDM')
ax.axhline(y=0,color='black',linewidth=0.5)
fig.set_size_inches(4,4)
plt.tight_layout()
plt.savefig('figures/HDDM_subj_v.png')

fig, ax = plt.subplots()
ax.bar(x_pos[4:11], means[4:11], yerr=errors[4:11], align='center', alpha=0.5, ecolor='black', capsize=10, width=0.38, color=['blue','red','blue','red','blue','red'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos[2:6])
ax.set_xticklabels(xtick_labels[2:6])
ax.set_title('Subj Rating informed HDDM')
ax.axhline(y=0,color='black',linewidth=0.5)
fig.set_size_inches(6,4)
plt.tight_layout()
plt.savefig('figures/HDDM_subj_zat.png')
