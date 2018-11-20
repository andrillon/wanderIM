#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 08:00:35 2018

@author: tand0009
"""
import pandas as pd
import matplotlib.pyplot as plt

data = hddm.load_csv('/Users/tand0009/anaconda2/envs/hddm/lib/python3.5/site-packages/hddm/examples/cavanagh_theta_nn.csv')

data = hddm.utils.flip_errors(data)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)

m_stim = hddm.HDDM(data, depends_on={'v': 'stim'})
m_stim.find_starting_values()
m_stim.sample(8000, burn=2000)

from patsy import dmatrix
dmatrix("C(stim, Treatment('WL'))", data.head(10))

m_within_subj = hddm.HDDMRegressor(data, "v ~ C(stim, Treatment('WL'))")

m_within_subj.sample(5000, burn=200)

v_WL, v_LL, v_WW = m_within_subj.nodes_db.ix[["v_Intercept",
                                              "v_C(stim, Treatment('WL'))[T.LL]",
                                              "v_C(stim, Treatment('WL'))[T.WW]"], 'node']
hddm.analyze.plot_posterior_nodes([v_WL, v_LL, v_WW])
