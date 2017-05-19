# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 12:41:05 2016

@author: aeber
"""
import scipy.stats as stats
import math

#returns an array of observation times where the ith element corresponds to
#detecting signal[i] given sum(noise[])[i] where both signal and noise are arrays
#of poisson distributed random variables, if passed key parameter pi then
#pi is an array of the signal to noise ratio (ie signal[i]/sum(noise[])[i])
def minimumTime(signal, *noise, **kwargs):
    T = []
    for i in range(len(signal)):
        S = signal[i]
        if "pi" in kwargs:
            pi = kwargs["pi"][i]
            T.append(25./(pi*S))
        else:
            N = 0
            for j in noise:
                N += noise[j][i]
                T.append(25.*float(N)/(S**2))
    return T


#returns a matrix of observation times where the ith element corresponds to 
#the observation time (m^2*s) required to distinguish the data[i] from model[i]
#assuming each element in data and model is a poisson distributed random variable
def distinguishTime(data, model):
    T = []
    for i in range(len(data)):
        d2 = (data[i] - model[i])**2
        T.append(50.*model[i]/d2)
    return T


#returns the observation time required to distinguish 
def distinguishTimeTotal(data, model):
    return distinguishTime([sum(data)], [sum(model)])[0]


#returns the probability that given data can be described by given model
#if using fluxes then T is the observation time
def chiSquared(data, model, T = 1):
    chi2 = 0
    df = 0
    for i in range(len(data)):
        if model[i] > 0:
            df += 1
            chi2 += (float(data[i] - model[i])**2)/model[i]
    chi2 *= T
    print chi2
    return 1. - stats.chi2.cdf(chi2, df)


#finds the -2ln(lambda) for a set of results n with expectations mu
#with fitted parameters m with observation time T
#if full is true then return the p value corresponding to these data being described by the same pdfs
def lnlikelyRatio(n, mu, m = 1, full = False):
    l, df = 0, 0
    for i in range(len(n)):
        if n[i] > 0 and mu[i] > 0:
            df += 1
            mu_ = mu[i]
            n_ = n[i]
            l += mu_ - n_ + n_*math.log(float(n_)/mu_)
    l *= 2.
    if full:
        return l, df, stats.chi2.cdf(l, df - m)
    else:
        return l
    

#poisson picks values for expectation value matrix mu    
def poissonPick(mu):
    events = []
    for i in range(len(mu)):
        events.append(stats.poisson.rvs(mu[i]))
    return events