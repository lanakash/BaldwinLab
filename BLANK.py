# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 10:31:14 2017

@author: Lana
"""

import scipy.stats as st

x = [0, 0.01]
y = [0, 0.02]

print st.mannwhitneyu(x,y)