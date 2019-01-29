# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 21:46:09 2017

@author: acburns

HDDM Script for DSPD Go/No-Go Data
"""
## Import required packages
import hddm
import numpy as np        
import pandas as pd
import matplotlib.pyplot as plt
from kabuki.analyze import gelman_rubin

# Load data from csv file into a NumPy structured array
data = mydata = hddm.load_csv('GNGData_270518_outrem.csv')

       
# Create histogram of RTs by subject, looks odd due to dummy coding of nogo response.
data = hddm.utils.flip_errors(data)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby(['session']):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)

plt.savefig('hddm_RThists_bysubj.pdf')


# HDDMStimCoding Model Subclass enables you to specify your bounds as go/no-go, as opposed to correct/incorrect. This allows the estimation of both z bias and v bias. 
# Code below fits the full model where session (condition) influences all parameters, tends to have the lowest DIC. 
model = hddm.HDDMStimCoding(data, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'cond_v', 'a': 'session', 't': 'session', 'z': 'session'}, p_outlier=.05)
model.find_starting_values()# Create model and start MCMC sampling
model.sample(30000, burn=10000, thin=2, dbname='hddm_stimZ&v_gonogo.db', db='pickle')
model.save('hddm_stimZ&v_gonogo')
model.print_stats()


# Plot posterior distributions and theoretical RT distributions
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

v_1, v_2 = model.nodes_db.node[['v(1)', 'v(2)']]
hddm.analyze.plot_posterior_nodes([v_1, v_2])
plt.xlabel('Drift')
plt.ylabel('Posterior probability')
plt.title('Posterior of Drift between timepoints')
plt.savefig('hddm_demo_fig_09.pdf')

print("P(v(go2) > v(go1)) = ", (v_1.trace() > v_2.trace()).mean())

dc_1, dc_2 = model.nodes_db.node[['dc(1)', 'dc(2)']]
hddm.analyze.plot_posterior_nodes([dc_1, dc_2])
plt.xlabel('dc')
plt.ylabel('Posterior probability')
plt.title('Posterior of dc between timepoints')
plt.savefig('hddm_demo_fig_09.pdf')

print("P(v(dc2) > v(dc1)) = ", (dc_1.trace() > dc_2.trace()).mean())

#Posterior predictive check
ppc_data = hddm.utils.post_pred_gen(model)
ppc_compare = hddm.utils.post_pred_stats(data, ppc_data)
ppc_stats = hddm.utils.post_pred_stats(data, ppc_data, call_compare=False)
print(ppc_compare)
ppc_data.to_csv(os.path.join(fig_dir, 'diagnostics', 'ppc_210518.csv'))

#Get Gelman-Rubin Stats
models = []
for i in range(5):
    m = hddm.HDDMStimCoding(data, include={'z'}, stim_col= 'stimulus', split_param='z', depends_on={'v': 'cond_v', 'a': 'session', 't': 'session', 'z': 'session'}, p_outlier=.05)
    m.find_starting_values()
    m.sample(5000, burn=2000)
    models.append(m)
model_base_name = '2018_GNG_HDDM'
fig_dir = os.path.join(model_base_name)

gelman = gelman_rubin(models)
gelman = pd.DataFrame.from_dict(gelman,orient='index')
gelman.to_csv(os.path.join(fig_dir, 'diagnostics', 'gelman_GNG.csv'))


##### Below code is saving data and creating figures

model_base_name = '2018_GNG_HDDM'

fig_dir = os.path.join(model_base_name)
results = model.gen_stats() # check the outputs in comparision with the true param
results.to_csv(os.path.join(fig_dir, 'diagnostics', 'results_ARC.csv'))

#data.to_csv(os.path.join(fig_dir, 'diagnostics', 'simdataTest.csv'))

text_file = open(os.path.join(fig_dir, 'diagnostics', 'DIC.txt'), 'w')
text_file.write("GNG_ARC".format(model, model.dic))
text_file.close()


params_of_interest = ['v(go1)', 'v(go2)', 'v(nogo1)', 'v(nogo2)', 't(1)', 't(2)', 'z(1)', 'z(2)', 'a(1)', 'a(2)']

traces = []
for p in range(len(params_of_interest)):
    traces.append(model.nodes_db.node[params_of_interest[p]].trace.gettrace())
    
# pass z values through link
traces[6] = 1-traces[6]
traces[7] = 1-traces[7] 


tracesarray = np.asarray(traces)
tracesFrame= pd.DataFrame(data=tracesarray) 
tracesFrame.to_csv(os.path.join(fig_dir, 'diagnostics', 'traces_ARC.csv'))

# Create stats
stats = []
for p in range(len(params_of_interest)):
    stat = np.min(np.mean(traces[p]))
    stats.append(stat)
stats = np.array(stats)
statsFrame= pd.DataFrame(data=stats[1:,])
statsFrame.to_csv(os.path.join(fig_dir, 'diagnostics', 'means_ARC.csv'))

