# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 19:01:47 2015

@author: cra
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# sanity check analytical solution

n0 = 3.
er = 2.
sigma = 4.

k = sqrt(n0/er)
kv = k*sqrt(2)


u0 = 1.0986122405839509

C = np.log(np.tanh(u0/4))

u_anl = linspace(0,u0,100)
x = u_anl.copy()
for p in range(np.alen(x)):
    x[p] = (C-np.log(np.tanh(u_anl[p]/4)))/kv

plt.plot(x,u_anl)
