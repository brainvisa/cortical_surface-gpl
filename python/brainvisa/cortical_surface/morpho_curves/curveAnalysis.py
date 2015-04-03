# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 15:14:51 2014

@author: olivier coulon
"""

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
def GetL1L2(signal, t):
    dt=0.05
    smooth=Diffusion(signal, dt, t)
    st=signal.std()
    minP=100
    maxP=-100
    l1=0
    l2=0
    for i in range(20, 50):
        if ((smooth[i]<minP) & (smooth[i]<=smooth[i+1]) & (smooth[i]<=smooth[i-1])):
            l1=i
            minP=smooth[i]
    for i in range(l1+1, 80):
        if ((smooth[i]>maxP) & (smooth[i]>=smooth[i+1]) & (smooth[i]>=smooth[i-1]) & ((smooth[i]-minP) > 0.5*st)):
            l2=i
            maxP=smooth[i]
            break
    return l1, l2, smooth
