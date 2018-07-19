#!/usr/bin/python

# Script by Wenyang Li

import numpy as np

a1 = np.array([15,15*np.sqrt(3)])
a2 = np.array([15,-15*np.sqrt(3)])
a3 = -a1-a2
d0 = np.array([20,0])
d1 = np.array([-10,10*np.sqrt(3)])
d2 = np.array([-10,-10*np.sqrt(3)])
n = 0
pos = {}
for ii in range(10):
    for jj in range(10):
        pos[n] = d0+ii*a1+jj*a2
        n += 1
for ii in range(11):
    for jj in range(10):
        pos[n] = d1+ii*a3+jj*a1
        n += 1
for ii in range(11):
    for jj in range(11):
        pos[n] = d2+ii*a2+jj*a3
        n += 1

print pos
print len(pos)

# Outriggers:

#for ii in range(4):
#    pos[n] = -30*a2 - ii*10*a3
#    n += 1
#for ii in range(5):
#    pos[n] = -20*a2 + 10*a3 - ii*10*a3
#    n += 1
#for ii in range(6):
#    if ii in [2,3]: continue
#    pos[n] = -10*a2 + 20*a3 - ii*10*a3
#    n += 1
#for ii in range(7):
#    if ii in [2,3,4]: continue
#    pos[n] = 30*a3 - ii*10*a3
#    n += 1
#for ii in range(6):
#    if ii in [2,3]: continue
#    pos[n] = 10*a2 + 30*a3 - ii*10*a3
#    n += 1
#for ii in range(5):
#    pos[n] = 20*a2 + 30*a3 - ii*10*a3
#    n += 1
#for ii in range(4):
#    pos[n] = 30*a2 +30*a3 - ii*10*a3
#    n += 1
#antpos = {'nant':n}
#for k in pos.keys():
#    antpos[k] = {}
#    antpos[k]['top_x']=pos[k][0]
#    antpos[k]['top_y']=pos[k][1]
