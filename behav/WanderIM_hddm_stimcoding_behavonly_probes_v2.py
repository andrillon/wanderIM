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
data = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults2_MW_new.txt')
data.columns = ['subj_idx','nblock','nprobe','task','ntrial','look','state','orig','awa','int','eng','perf','vig','corr','rt','stim','dprobe','response','stimulus','vigC','cond_v']
#data = hddm.load_csv('/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults2.txt')
#data.columns = ['subj_idx','nblock','task','ntrial','stimid','correctness','rt','stimulus','response','cond_v']
data = data[np.logical_or(np.logical_and(data.dprobe > -19,data.stimulus ==1),np.logical_and(data.dprobe > -3,data.stimulus ==0))]
data.fillna(999, inplace=True)
#data = data[data.rt > 0.2]

# Create histogram of RTs by subject, looks odd due to dummy coding of nogo response.
data = hddm.utils.flip_errors(data)
data = data[data.state != 4]
data["state"] = data["state"].astype('category')

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby(['state']):
    subj_data.rt.hist(bins=50, histtype='step', ax=ax)

#plt.savefig('hddm_RThists_bysubj.pdf')

dataD = data[data.task == 2]
dataF = data[data.task == 1]

# HDDMStimCoding Model Subclass enables you to specify your bounds as go/no-go, as opposed to correct/incorrect. This allows the estimation of both z bias and v bias. 
# Code below fits the full model where session (condition) influences all parameters, tends to have the lowest DIC. 
modelD = hddm.HDDMStimCoding(dataD, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'cond_v', 'a': 'state', 'z': 'state', 't': 'state'}, p_outlier=.05)
modelD.find_starting_values()# Create model and start MCMC sampling
modelD.sample(5000, burn=200, dbname='hddm_stim.db', db='pickle')
#model.save('hddm_stim')
modelD.print_stats()

modelD.plot_posteriors(save=False)
modelD.plot_posterior_predictive()
modelD.plot_posteriors_conditions()

# HDDMStimCoding Model Subclass enables you to specify your bounds as go/no-go, as opposed to correct/incorrect. This allows the estimation of both z bias and v bias. 
# Code below fits the full model where session (condition) influences all parameters, tends to have the lowest DIC. 
modelF = hddm.HDDMStimCoding(dataF, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'cond_v', 'a': 'state', 'z': 'state', 't': 'state'}, p_outlier=.05)
modelF.find_starting_values()# Create model and start MCMC sampling
modelF.sample(5000, burn=200, dbname='hddm_stim.db', db='pickle')
#model.save('hddm_stim')
modelF.print_stats()

modelF.plot_posteriors(save=False)
modelF.plot_posterior_predictive()
modelF.plot_posteriors_conditions()

# Plot posterior probabilities for parameters and determine overlap
z_1, z_2, z_3 = modelD.nodes_db.node[['z(1)', 'z(2)', 'z(3)']]
hddm.analyze.plot_posterior_nodes([z_1, z_2, z_3])
plt.xlabel('Response Bias z')
plt.ylabel('Posterior probability')
plt.title('Posterior of response bias between timepoints')

print("P(z(2) > z(1)) = ", (z_1.trace() > z_2.trace()).mean())
print("P(z(3) > z(1)) = ", (z_1.trace() > z_3.trace()).mean())
print("P(z(3) > z(2)) = ", (z_2.trace() > z_3.trace()).mean())

a_1, a_2, a_3 = modelD.nodes_db.node[['a(1)', 'a(2)', 'a(3)']]
hddm.analyze.plot_posterior_nodes([a_1, a_2, a_3])
plt.xlabel('Boundary setting a')
plt.ylabel('Posterior probability')
plt.title('Posterior of boundary setting between timepoints')

print("P(a(2) > a(1)) = ", (a_1.trace() > a_2.trace()).mean())
print("P(a(3) > a(1)) = ", (a_1.trace() > a_3.trace()).mean())
print("P(a(3) > a(2)) = ", (a_2.trace() > a_3.trace()).mean())

t_1, t_2, t_3 = modelD.nodes_db.node[['t(1)', 't(2)', 't(3)']]
hddm.analyze.plot_posterior_nodes([t_1, t_2, t_3])
plt.xlabel('NDT')
plt.ylabel('Posterior probability')
plt.title('Posterior of NDT between timepoints')

print("P(t(2) < t(1)) = ", (t_1.trace() < t_2.trace()).mean())
print("P(t(3) > t(1)) = ", (t_1.trace() > t_3.trace()).mean())
print("P(t(3) > t(2)) = ", (t_2.trace() > t_3.trace()).mean())


v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3 = modelD.nodes_db.node[['v(go_on)', 'v(go_mw)', 'v(go_mb)','v(nogo_on)', 'v(nogo_mw)', 'v(nogo_mb)']]
hddm.analyze.plot_posterior_nodes([v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3])
plt.xlabel('Drift')
plt.ylabel('Posterior probability')
plt.title('Posterior of Drift between timepoints')

print("P(v(go_on) > v(go_mw)) = ", (v_g1.trace() > v_g2.trace()).mean())
print("P(v(go_on) > v(go_mb)) = ", (v_g1.trace() > v_g3.trace()).mean())
print("P(v(go_mw) > v(go_mb)) = ", (v_g2.trace() > v_g3.trace()).mean())
print("P(v(nogo_on) > v(nogo_mw)) = ", (v_ng1.trace() > v_ng2.trace()).mean())
print("P(v(nogo_on) > v(nogo_mb)) = ", (v_ng1.trace() > v_ng3.trace()).mean())
print("P(v(nogo_mw) > v(nogo_mb)) = ", (v_ng2.trace() > v_ng3.trace()).mean())

v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3 = modelD.nodes_db.node[['v(go_on)', 'v(go_mw)', 'v(go_mb)','v(nogo_on)', 'v(nogo_mw)', 'v(nogo_mb)']]
#t_1, t_2, t_3 = modelD.nodes_db.node[['t(1)', 't(2)', 't(3)']]
a_1, a_2, a_3 = modelD.nodes_db.node[['a(1)', 'a(2)', 'a(3)']]
z_1, z_2, z_3 = modelD.nodes_db.node[['z(1)', 'z(2)', 'z(3)']]
mean_v_g1=(v_g1.trace()).mean()
mean_v_g2=(v_g2.trace()).mean()
mean_v_g3=(v_g3.trace()).mean()
mean_v_ng1=(v_ng1.trace()).mean()
mean_v_ng2=(v_ng2.trace()).mean()
mean_v_ng3=(v_ng3.trace()).mean()
mean_z_1=(z_1.trace()).mean()
mean_z_2=(z_2.trace()).mean()
mean_z_3=(z_3.trace()).mean()
mean_a_1=(a_1.trace()).mean()
mean_a_2=(a_2.trace()).mean()
mean_a_3=(a_3.trace()).mean()
#mean_t_1=(t_1.trace()).mean()
#mean_t_2=(t_2.trace()).mean()
#mean_t_3=(t_3.trace()).mean()

std_v_g1=(v_g1.trace()).std()
std_v_g2=(v_g2.trace()).std()
std_v_g3=(v_g3.trace()).std()
std_v_ng1=(v_ng1.trace()).std()
std_v_ng2=(v_ng2.trace()).std()
std_v_ng3=(v_ng3.trace()).std()
std_z_1=(z_1.trace()).std()
std_z_2=(z_2.trace()).std()
std_z_3=(z_3.trace()).std()
std_a_1=(a_1.trace()).std()
std_a_2=(a_2.trace()).std()
std_a_3=(a_3.trace()).std()
#std_t_1=(t_1.trace()).std()
#std_t_2=(t_2.trace()).std()
#std_t_3=(t_3.trace()).std()

x_pos=np.array([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5])+np.array([-0.2,0,0.2,-0.2,0,0.2,-0.2,0,0.2,-0.2,0,0.2,-0.2,0,0.2])
xtick_pos=np.array([1,2,3,4,5])
xtick_labels = ['v_{go}', 'v_{nogo}', 'z', 'a', 't']
meansD = [mean_v_g1, mean_v_g2, mean_v_g3, mean_v_ng1, mean_v_ng2, mean_v_ng3, mean_a_1, mean_a_2, mean_a_3, mean_z_1, mean_z_2, mean_z_3]#, mean_t_1, mean_t_2, mean_t_3]
errorsD = [std_v_g1, std_v_g2, std_v_g3, std_v_ng1, std_v_ng2, std_v_ng3, std_a_1, std_a_2, std_a_3, std_z_1, std_z_2, std_z_3]#, std_t_1, std_t_2, std_t_3]

#fig, ax = plt.subplots()
#ax.bar(x_pos, means, yerr=errors, align='center', alpha=0.5, ecolor='black', capsize=10, width=0.18, color=['green','orange','blue','green','orange','blue','green','orange','blue','green','orange','blue','green','orange','blue'])
#ax.set_ylabel('HDDM Parameters')
#ax.set_xticks(xtick_pos)
#ax.set_xticklabels(xtick_labels)
#ax.axhline(y=0,color='black',linewidth=0.5)
#ax.set_title('Local Sleep informed HDDM')

fig, ax = plt.subplots()
ax.bar(x_pos[0:9], meansD[0:9], yerr=errorsD[0:9], align='center', alpha=0.5, ecolor='black', capsize=10, width=0.18, color=['green','orange','blue','green','orange','blue','green','orange','blue'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos[0:3])
ax.set_xticklabels(xtick_labels[0:3])
ax.set_title('Local Sleep informed HDDM')
ax.axhline(y=0,color='black',linewidth=0.5)
plt.tight_layout()

fig, ax = plt.subplots()
ax.bar(x_pos[9:12], meansD[9:12], yerr=errorsD[9:12], align='center', alpha=0.5, ecolor='black', capsize=10, width=0.18, color=['green','orange','blue','green','orange','blue'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos[3:4])
ax.set_xticklabels(xtick_labels[3:4])
ax.set_title('Local Sleep informed HDDM')
ax.axhline(y=0,color='black',linewidth=0.5)
plt.tight_layout()

v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3 = modelD.nodes_db.node[['v(go_on)', 'v(go_mw)', 'v(go_mb)','v(nogo_on)', 'v(nogo_mw)', 'v(nogo_mb)']]
fig = plt.figure()
dc1=(v_g1.trace() + v_ng1.trace())
dc2=(v_g2.trace() + v_ng2.trace())
dc3=(v_g3.trace() + v_ng3.trace())
plt.hist(dc1,color='green',bins=50,alpha=0.5);
plt.hist(dc2,color='orange',bins=50,alpha=0.5);
plt.hist(dc3,color='blue',bins=50,alpha=0.5);

# Plot posterior probabilities for parameters and determine overlap
z_1, z_2, z_3 = modelF.nodes_db.node[['z(1)', 'z(2)', 'z(3)']]
hddm.analyze.plot_posterior_nodes([z_1, z_2, z_3])
plt.xlabel('Response Bias z')
plt.ylabel('Posterior probability')
plt.title('Posterior of response bias between timepoints')

print("P(z(2) > z(1)) = ", (z_1.trace() > z_2.trace()).mean())
print("P(z(3) > z(1)) = ", (z_1.trace() > z_3.trace()).mean())
print("P(z(3) > z(2)) = ", (z_2.trace() > z_3.trace()).mean())

a_1, a_2, a_3 = modelF.nodes_db.node[['a(1)', 'a(2)', 'a(3)']]
hddm.analyze.plot_posterior_nodes([a_1, a_2, a_3])
plt.xlabel('Boundary setting a')
plt.ylabel('Posterior probability')
plt.title('Posterior of boundary setting between timepoints')

print("P(a(2) > a(1)) = ", (a_1.trace() > a_2.trace()).mean())
print("P(a(3) > a(1)) = ", (a_1.trace() > a_3.trace()).mean())
print("P(a(3) > a(2)) = ", (a_2.trace() > a_3.trace()).mean())

t_1, t_2, t_3 = modelF.nodes_db.node[['t(1)', 't(2)', 't(3)']]
hddm.analyze.plot_posterior_nodes([t_1, t_2, t_3])
plt.xlabel('NDT')
plt.ylabel('Posterior probability')
plt.title('Posterior of NDT between timepoints')

print("P(t(2) < t(1)) = ", (t_1.trace() < t_2.trace()).mean())
print("P(t(3) > t(1)) = ", (t_1.trace() > t_3.trace()).mean())
print("P(t(3) > t(2)) = ", (t_2.trace() > t_3.trace()).mean())


v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3 = modelF.nodes_db.node[['v(go_on)', 'v(go_mw)', 'v(go_mb)','v(nogo_on)', 'v(nogo_mw)', 'v(nogo_mb)']]
hddm.analyze.plot_posterior_nodes([v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3])
plt.xlabel('Drift')
plt.ylabel('Posterior probability')
plt.title('Posterior of Drift between timepoints')

print("P(v(go_on) > v(go_mw)) = ", (v_g1.trace() > v_g2.trace()).mean())
print("P(v(go_on) > v(go_mb)) = ", (v_g1.trace() > v_g3.trace()).mean())
print("P(v(go_mw) > v(go_mb)) = ", (v_g2.trace() > v_g3.trace()).mean())
print("P(v(nogo_on) > v(nogo_mw)) = ", (v_ng1.trace() > v_ng2.trace()).mean())
print("P(v(nogo_on) > v(nogo_mb)) = ", (v_ng1.trace() > v_ng3.trace()).mean())
print("P(v(nogo_mw) > v(nogo_mb)) = ", (v_ng2.trace() > v_ng3.trace()).mean())

v_g1, v_g2, v_g3, v_ng1, v_ng2, v_ng3 = modelF.nodes_db.node[['v(go_on)', 'v(go_mw)', 'v(go_mb)','v(nogo_on)', 'v(nogo_mw)', 'v(nogo_mb)']]
#t_1, t_2, t_3 = modelF.nodes_db.node[['t(1)', 't(2)', 't(3)']]
a_1, a_2, a_3 = modelF.nodes_db.node[['a(1)', 'a(2)', 'a(3)']]
z_1, z_2, z_3 = modelF.nodes_db.node[['z(1)', 'z(2)', 'z(3)']]
mean_v_g1=(v_g1.trace()).mean()
mean_v_g2=(v_g2.trace()).mean()
mean_v_g3=(v_g3.trace()).mean()
mean_v_ng1=(v_ng1.trace()).mean()
mean_v_ng2=(v_ng2.trace()).mean()
mean_v_ng3=(v_ng3.trace()).mean()
mean_z_1=(z_1.trace()).mean()
mean_z_2=(z_2.trace()).mean()
mean_z_3=(z_3.trace()).mean()
mean_a_1=(a_1.trace()).mean()
mean_a_2=(a_2.trace()).mean()
mean_a_3=(a_3.trace()).mean()
#mean_t_1=(t_1.trace()).mean()
#mean_t_2=(t_2.trace()).mean()
#mean_t_3=(t_3.trace()).mean()

std_v_g1=(v_g1.trace()).std()
std_v_g2=(v_g2.trace()).std()
std_v_g3=(v_g3.trace()).std()
std_v_ng1=(v_ng1.trace()).std()
std_v_ng2=(v_ng2.trace()).std()
std_v_ng3=(v_ng3.trace()).std()
std_z_1=(z_1.trace()).std()
std_z_2=(z_2.trace()).std()
std_z_3=(z_3.trace()).std()
std_a_1=(a_1.trace()).std()
std_a_2=(a_2.trace()).std()
std_a_3=(a_3.trace()).std()
#std_t_1=(t_1.trace()).std()
#std_t_2=(t_2.trace()).std()
#std_t_3=(t_3.trace()).std()

x_pos=np.array([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5])+np.array([-0.2,0,0.2,-0.2,0,0.2,-0.2,0,0.2,-0.2,0,0.2,-0.2,0,0.2])
xtick_pos=np.array([1,2,3,4,5])
xtick_labels = ['v_{go}', 'v_{nogo}', 'z', 'a', 't']
meansD = [mean_v_g1, mean_v_g2, mean_v_g3, mean_v_ng1, mean_v_ng2, mean_v_ng3, mean_a_1, mean_a_2, mean_a_3, mean_z_1, mean_z_2, mean_z_3]#, mean_t_1, mean_t_2, mean_t_3]
errorsD = [std_v_g1, std_v_g2, std_v_g3, std_v_ng1, std_v_ng2, std_v_ng3, std_a_1, std_a_2, std_a_3, std_z_1, std_z_2, std_z_3]#, std_t_1, std_t_2, std_t_3]

#fig, ax = plt.subplots()
#ax.bar(x_pos, means, yerr=errors, align='center', alpha=0.5, ecolor='black', capsize=10, width=0.18, color=['green','orange','blue','green','orange','blue','green','orange','blue','green','orange','blue','green','orange','blue'])
#ax.set_ylabel('HDDM Parameters')
#ax.set_xticks(xtick_pos)
#ax.set_xticklabels(xtick_labels)
#ax.axhline(y=0,color='black',linewidth=0.5)
#ax.set_title('Local Sleep informed HDDM')

fig, ax = plt.subplots()
ax.bar(x_pos[0:9], meansD[0:9], yerr=errorsD[0:9], align='center', alpha=0.5, ecolor='black', capsize=10, width=0.18, color=['green','orange','blue','green','orange','blue','green','orange','blue'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos[0:3])
ax.set_xticklabels(xtick_labels[0:3])
ax.set_title('Local Sleep informed HDDM')
ax.axhline(y=0,color='black',linewidth=0.5)
plt.tight_layout()

fig, ax = plt.subplots()
ax.bar(x_pos[9:12], meansD[9:12], yerr=errorsD[9:12], align='center', alpha=0.5, ecolor='black', capsize=10, width=0.18, color=['green','orange','blue','green','orange','blue'])
ax.set_ylabel('HDDM Parameters')
ax.set_xticks(xtick_pos[3:4])
ax.set_xticklabels(xtick_labels[3:4])
ax.set_title('Local Sleep informed HDDM')
ax.axhline(y=0,color='black',linewidth=0.5)
plt.tight_layout()

fig = plt.figure()
dc1=(v_g1.trace() + v_ng1.trace())
dc2=(v_g2.trace() + v_ng2.trace())
dc3=(v_g3.trace() + v_ng3.trace())
plt.hist(dc1,color='green',bins=50,alpha=0.5);
plt.hist(dc2,color='orange',bins=50,alpha=0.5);
plt.hist(dc3,color='blue',bins=50,alpha=0.5);
