#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:22:08 2019

@author: simon
"""

import numpy as np
import numpy.ma as ma
import copy, os, sys
from PIL import Image
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.path as mlp
from scipy.optimize import curve_fit


import FPF_PeakFourierFunctions as PFF




print '# ============================================================================='
print '# Test cases for components of PeakFourierFits'
print '# =============================================================================\n'

print 'The set of d0s, Hs, Ws, Ps (for three peaks) and a quartic backgound is used for the test cases:'


d0s = [[3.34, .12321, 0.2], [4.5], [6.1234]]
Hs = [[12.3, 2.005, -0.04], [8.04], [10.2]]
Ws = [[.1], [.05, .1, .1], [.4]]
Ps = [[1], [0], [0.534, .05, -.05]]

bg = [[4,.3, .3], [-.9], [.3], [.05, .01, .01]]

print '   d0s are:', d0s
print '   Hs are: ', Hs
print '   Ws are: ', Ws
print '   Ps are: ', Ps
print '   background is:', bg


# =============================================================================
# Test GuessGet and GuessApply
# =============================================================================
print '\nTest 1: comfirm GuessGet and GuessApply work.'
print '  These functions apply a ''ChangeArray'' to the arrays of parameters to get values for fitting'
print '  and paste fitted value back into the arrays.'

print '  GuessGet extracts the values associated with the non-zero values in the ''ChangeArray''.'

ChangeArr = [[[0,3,0,0],[0,0,0,0],[0,0,0,0]],[0,0,0,0]]
print '    This ''ChangeArray'' acts on the heights of the first peak.'
print '    ', ChangeArr
p0 = PFF.GuessGet(ChangeArr, d0s, Hs, Ws, Ps, bg)
print '    The values it returns are:'
print '        ', p0
print '    Which should be the same as ''Hs[0]'':'''
print '        ', Hs[0]
if not p0 == Hs[0]:
    sys.exit('Failed test: the arrays are not the same')
else:
    print '      PASSED\n'

ChangeArr = [[[0,0,0,0],[0,0,0,0],[0,0,0,0]],[3,1,1,3]]
print '    This ''ChangeArray'' gets all the background parameters as an array.'
print '    ', ChangeArr
p0 = PFF.GuessGet(ChangeArr, d0s, Hs, Ws, Ps, bg)
print '    The values it returns are:'
print '        ', p0
print '    Which should be the same as the background values in a single array:'
print '        ', [item for sublist in bg for item in sublist]
if not p0 == [item for sublist in bg for item in sublist]:
    sys.exit('Failed test: the arrays are not the same')
else:
    print '      PASSED\n'    
    

print '  GuessApply takes an array of fitted values and the same ''changearray'' and pastes the values'
print '  into the paramter arrays.'  

ChangeArr = [[[0,3,0,0],[0,0,0,0],[0,0,0,0]],[0,0,0,0]]
print '    Using ''ChangeArray'' that acts on the heights of the first peak.'
print '    ', ChangeArr
print '    We replace the values in Hs with new values:'
preturn = [6.,7.,8.]
print '    ', preturn
print '    Hs before:'
print '    ', Hs
d0s_new, Hs_new, Ws_new, Ps_new, bg_new = PFF.GuessApply(ChangeArr, d0s, Hs, Ws, Ps, bg, preturn)
print '    New Hs are:'
print '    ', Hs_new
print '    d0s are:'
print '    ', d0s_new
print '    Ws are:'
print '    ', Ws_new
print '    Ps are:'
print '    ', Ps_new
print '    background is:'
print '    ', bg_new
if not (d0s==d0s_new and Ws==Ws_new and Ps==Ps_new and bg==bg_new and Hs_new[0]==preturn and Hs[1:]==Hs_new[1:]):
    sys.exit('Failed test: the arrays are not the same')
else:
    print '      PASSED\n'

ChangeArr = [[[0,0,0,0],[0,0,0,0],[0,0,0,0]],[3,1,1,3]]
print '    Using ''ChangeArray'' that acts on the background.'
print '    ', ChangeArr
print '    We replace the values in bg with new values:'
preturn = [-4,-.3, -.3, .9, -.3, -.05, -.01, -.01]
print '    ', preturn
print '    background before:'
bg_in =  bg[:]
print '    ', bg_in
d0s_new, Hs_new, Ws_new, Ps_new, bg_new = PFF.GuessApply(ChangeArr, d0s, Hs, Ws, Ps, bg_in, preturn)
print 'Oh dear', bg
print '    New bg is:'
print '    ', bg_new
print '    d0s are:'
print '    ', d0s_new
print '    Hs are:'
print '    ', Hs_new
print '    Ws are:'
print '    ', Ws_new
print '    Ps are:'
print '    ', Ps_new
if not (d0s==d0s_new and Ws==Ws_new and Ps==Ps_new and Hs==Hs_new and [item for sublist in bg_new for item in sublist]==preturn):
    sys.exit('Failed test: the arrays are not the same')
else:
    print '      PASSED\n'


# =============================================================================
# Single Peak on quartic background
# =============================================================================
print '\nTest 2: ''Test PeaksModel'' - Make 1D, Single Peak on quartic background and plot.'

twotheta = np.linspace(2,5,400)
azi = twotheta * 0
print '  Range of plot:', np.min(twotheta), 'to',  np.max(twotheta)
d0s_now = [d0s[0][0]]
print '  Peak centered on', d0s_now
Hs_now  = [Hs[0][0]]
print '  Peak height     ', Hs_now
Ws_now  = [Ws[0][0]]
print '  Peak width      ', Ws_now
Ps_now  = [Ps[0][0]]
print '  Peak profile    ', Ps_now
bg_now  = [[bg[0][0]],[bg[1][0]],[bg[2][0]],[bg[3][0]]]
print '  background      ', bg_now


mod, peaks = PFF.PeaksModel(twotheta, azi, d0s_now, Hs_now, Ws_now, Ps_now, bg_now)

print 'Size input array: ', np.size(twotheta)
print 'Size output array:', np.size(mod)

plt.plot(twotheta, mod,'.-')
plt.show()

        
    
# =============================================================================
# 3 Peaks on quartic background
# =============================================================================
print '\nTest 3: ''Test PeaksModel'' - Make 1D, 3 peaks Single Peak on quartic background and plot.'

twotheta = np.linspace(2,8,400)
azi = twotheta * 0
print '  Range of plot:', np.min(twotheta), 'to',  np.max(twotheta)
d0s_now = [[d0s[0][0]],[d0s[1][0]],[d0s[2][0]]]
print '  Peak centered on', d0s_now
Hs_now  = [[Hs[0][0]],[Hs[1][0]],[Hs[2][0]]]
print '  Peak height     ', Hs_now
Ws_now  = [[Ws[0][0]],[Ws[1][0]],[Ws[2][0]]]
print '  Peak width      ', Ws_now
Ps_now  = [[Ps[0][0]],[Ps[1][0]],[Ps[2][0]]]
print '  Peak profile    ', Ps_now
bg_now  = [[bg[0][0]],[bg[1][0]],[bg[2][0]],[bg[3][0]]]
print '  background      ', bg_now

mod, peaks = PFF.PeaksModel(twotheta, azi, d0s_now, Hs_now, Ws_now, Ps_now, bg_now)

print 'Size input array:', np.size(twotheta)
print 'Size output array:', np.size(mod)

plt.plot(twotheta, mod,'.-')
plt.show()


   
# =============================================================================
# 2D test -- 3 Peaks on quartic background
# =============================================================================
print '\nTest 4: ''Test PeaksModel'' - Make 2D, 3 peaks Single Peak on quartic background and plot.'


twotheta = np.linspace(2,8,400)
azi      = np.linspace(0,360,360)
twotheta, azi   = np.meshgrid(twotheta, azi)

print '  Range of plot:', np.min(twotheta), 'to',  np.max(twotheta)
d0s_now = d0s
print '  Peak centered on', d0s_now
Hs_now  = Hs
print '  Peak height     ', Hs_now
Ws_now  = Ws
print '  Peak width      ', Ws_now
Ps_now  = Ps
print '  Peak profile    ', Ps_now
bg_now  = bg
print '  background      ', bg_now



mod, peaks = PFF.PeaksModel(twotheta, azi, d0s, Hs, Ws, Ps, bg)

print 'Size input array:', np.size(twotheta)
print 'Size output array:', np.size(mod)

plt.scatter(twotheta, azi, s=2, c=mod, edgecolors='none')
plt.show()





# =============================================================================
# 1D test -- Fit single peak -- all parameters
# =============================================================================
print '\nTest 5: ''Test FitModel'' - Fit single peak to linear background and plot.'


twotheta = np.linspace(2,5,400)
azi = twotheta * 0
print '  Range of plot:', np.min(twotheta), 'to',  np.max(twotheta)
d0s_now = [[d0s[0][0]]]
print '  Peak centered on', d0s_now
Hs_now  = [[Hs[0][0]]]
print '  Peak height     ', Hs_now
Ws_now  = [[Ws[0][0]]]
print '  Peak width      ', Ws_now
Ps_now  = [[Ps[0][0]]]
print '  Peak profile    ', Ps_now
bg_now  = [[bg[0][0]],[bg[1][0]]]
print '  background      ', bg_now
Conv    = np.array([0])

mod, peaks = PFF.PeaksModel(twotheta, azi, d0s_now, Hs_now, Ws_now, Ps_now, bg_now)

mod = np.random.normal(loc=0.0, scale=1.5, size=mod.size) + mod
ChangeArray = [[[1,1,1,1]],[1,1]]

print d0s_now
d0s_guess = np.array(d0s_now) + np.random.normal(loc=0.0, scale=.05, size=1)
print d0s_guess

d0s_new, Hs_new, Ws_new, Ps_new, bg_new,a,b = PFF.FitModel(mod, twotheta, azi, ChangeArray, d0s_guess[:], Hs_now[:], Ws_now[:], Ps_now[:], bg_now[:])


modfit, peaksfit = PFF.PeaksModel(twotheta, azi, d0s_new, Hs_new, Ws_new, Ps_new, bg_new)

print '  Fitted Params'
d0s_now = [[d0s[0][0]]]
print '  Peak centered on', d0s_new
print '  Peak height     ', Hs_new
print '  Peak width      ', Ws_new
print '  Peak profile    ', Ps_new
print '  background      ', bg_new


print 'Size input array: ', np.size(twotheta)
print 'Size output array:', np.size(mod)

plt.plot(twotheta, mod,'r.', twotheta, modfit,'b-',)
plt.show()






# =============================================================================
# 1D test -- Fit 2 peaks -- all parameters
# =============================================================================
print '\nTest 6: ''Test FitModel'' - Fit 2 peaks to quardractic background and plot.'


twotheta = np.linspace(2,5,400)
azi = twotheta * 0
print '  Range of plot:', np.min(twotheta), 'to',  np.max(twotheta)
d0s_now = [[d0s[0][0]], [d0s[0][0]+0.3]]
print '  Peak centered on', d0s_now
Hs_now  = [[Hs[0][0]], [Hs[1][0]]]
print '  Peak height     ', Hs_now
Ws_now  = [[Ws[0][0]],[Ws[1][0]]]
print '  Peak width      ', Ws_now
Ps_now  = [[Ps[0][0]],[Ps[1][0]]]
print '  Peak profile    ', Ps_now
bg_now  = [[bg[0][0]],[bg[1][0]],[bg[2][0]]]
print '  background      ', bg_now
Conv    = np.array([0])

mod, peaks = PFF.PeaksModel(twotheta, azi, d0s_now, Hs_now, Ws_now, Ps_now, bg_now)

mod = np.random.normal(loc=0.0, scale=1.5, size=mod.size) + mod
ChangeArray = [[[1,1,1,1], [1,1,1,1]],[1,1,1]]

print d0s_now
d0s_guess = np.array(d0s_now) + np.random.normal(loc=0.0, scale=.05, size=1)
print d0s_guess

d0s_new, Hs_new, Ws_new, Ps_new, bg_new,a,b = PFF.FitModel(mod, twotheta, azi, ChangeArray, d0s_guess[:], Hs_now[:], Ws_now[:], Ps_now[:], bg_now[:])


modfit, peaksfit = PFF.PeaksModel(twotheta, azi, d0s_new, Hs_new, Ws_new, Ps_new, bg_new)

print '  Fitted Params'
d0s_now = [[d0s[0][0]]]
print '  Peak centered on', d0s_new
print '  Peak height     ', Hs_new
print '  Peak width      ', Ws_new
print '  Peak profile    ', Ps_new
print '  background      ', bg_new


print 'Size input array: ', np.size(twotheta)
print 'Size output array:', np.size(mod)

plt.plot(twotheta, mod,'r.', twotheta, modfit,'b-',)
plt.show()







# =============================================================================
# 1D test -- Fit 2 peaks -- all parameters
# =============================================================================
print '\nTest 7: ''Test FitModel'' - Fit 2 peaks to 2D background and plot.'


twotheta = np.linspace(2,5,400)
azi      = np.linspace(0,360,360)
twotheta, azi   = np.meshgrid(twotheta, azi)

print '  Range of plot:', np.min(twotheta), 'to',  np.max(twotheta)
d0s_now = [d0s[0][0:3], [d0s[0][0]+0.3]]
print '  Peak centered on', d0s_now
Hs_now  = [[Hs[0][0]], [Hs[1][0]]]
print '  Peak height     ', Hs_now
Ws_now  = [[Ws[0][0]],[Ws[1][0]]]
print '  Peak width      ', Ws_now
Ps_now  = [[Ps[0][0]],[Ps[1][0]]]
print '  Peak profile    ', Ps_now
bg_now  = [[bg[0][0]],[bg[1][0]],[bg[2][0]]]
print '  background      ', bg_now
Conv    = np.array([0])

mod, peaks = PFF.PeaksModel(twotheta, azi, d0s_now, Hs_now, Ws_now, Ps_now, bg_now)
mod = np.random.normal(loc=0.0, scale=2, size=mod.shape) + mod
ChangeArray = [[[3,1,1,1], [1,1,1,1]],[1,1,1]]

print d0s_now
d0s_guess = [[3.3,0.1,0.1],[3.6]] #np.array(d0s_now) #+ np.random.normal(loc=0.0, scale=.05, size=1)
print d0s_guess

d0s_new, Hs_new, Ws_new, Ps_new, bg_new,a,b = PFF.FitModel(mod.flatten(), twotheta.flatten(), azi.flatten(), ChangeArray, d0s_guess[:], Hs_now[:], Ws_now[:], Ps_now[:], bg_now[:])


modfit, peaksfit = PFF.PeaksModel(twotheta, azi, d0s_new, Hs_new, Ws_new, Ps_new, bg_new)

print '  Fitted Params'
d0s_now = [[d0s[0][0]]]
print '  Peak centered on', d0s_new
print '  Peak height     ', Hs_new
print '  Peak width      ', Ws_new
print '  Peak profile    ', Ps_new
print '  background      ', bg_new


print 'Size input array: ', np.size(twotheta)
print 'Size output array:', np.size(mod)

fig = plt.figure()
ax = fig.add_subplot(1,2,1)
ax1 = plt.subplot(121)
plt.scatter(twotheta, azi, s=2, c=mod, edgecolors='none')
ax1.set_title('Data')
ax2 = plt.subplot(122)
plt.scatter(twotheta, azi, s=2, c=modfit, edgecolors='none')
ax2.set_title('Fit')
plt.show()
























