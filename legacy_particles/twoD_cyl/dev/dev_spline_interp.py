# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 12:04:07 2016

@author: cra

Test constrained spline representation of sparse function

SO: 32046582

"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')
import numpy as np


from scipy.interpolate import UnivariateSpline, splev, splrep


n = 5
#x = np.linspace(0, 2*np.pi, n)
x = np.asarray([0,0.5,1,1.5,2,4,10])

y0 = np.cos(x) # zero initial slope
std = 0.5
noise = np.random.normal(0, std, len(x))
y = y0 #+ noise

sp = fh.spline_con_gradient(x, y, [10],[0])

plt.figure()
X = np.linspace(x.min(), x.max(), len(x)*10)
plt.plot(X, sp(X), '-r', lw=2, label='spline')
plt.plot(X, sp.derivative()(X), '-g', label='slope')
plt.plot(x, y, 'ok', label='data')
plt.legend(loc='best')
plt.show()

fh.place_current_fig()