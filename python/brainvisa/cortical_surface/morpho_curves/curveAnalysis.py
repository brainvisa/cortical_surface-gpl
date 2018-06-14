# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 15:14:51 2014

@author: olivier coulon
"""
from __future__ import print_function

from numpy import *

# Laplacian of a 1D signal
def Laplacien(signal):
    s=signal.size
    lap=ones(s)
    lap[0]=0
    for i in range(1, s-1):
        lap[i]=signal[i-1]-2*signal[i]+signal[i+1]
    lap[s-1]=0
    return lap

# 1-D smoothing using Laplacian
def Diffusion(signal, dt, t):
    diff=signal
    n=int(t/dt)
    for i in linspace(0,t, n, endpoint=False):
        diff=diff + dt*0.5*Laplacien(diff)
    return diff

# When the curve is a sulcus profile, get L1 L2 
def GetL1L2(signal, t, species='human'):
    if species=='human':
        boundInf=50
        boundSup=90
    elif species=='baboon':
        boundInf=60
        boundSup=90
    else:
        boundInf=50
        boundSup=90
    dt=0.05
    smooth=Diffusion(signal, dt, t)
    st=signal.std()

    print('signalSTD=', st)
    print('signal range=', signal.max() - signal.min())
    minP=100
    maxP=-100
    l1=0
    l2=0
    for i in range(20, boundInf):
        if ((signal[i]<minP) & (smooth[i]<=smooth[i+1]) & (smooth[i]<=smooth[i-1])):
            l1=i
            minP=signal[i]
    for i in range(l1+1, boundSup):
#        if ((signal[i]>maxP) & (smooth[i]>=smooth[i+1]) & (smooth[i]>=smooth[i-1]) & ((signal[i]-minP) > 0.1*st)):
        if ( (smooth[i]>=smooth[i+1]) & (smooth[i]>=smooth[i-1])) :# & ((signal[i]-minP) > 0.1*st)):

#
            # nextmin=99
            # for j in range(i+1, boundSup):
            #      if ((smooth[j]<=smooth[j+1]) & (smooth[j]<=smooth[j-1])):
            #         nextmin=signal[j]
            #         j=boundSup
#            if ((signal[i]-nextmin) > 0.1*st):
#
            l2=i
            print('l2=', l2)
            maxP=smooth[i]
            i=boundSup

    print('l2=', l2)

    return l1, l2, smooth
 