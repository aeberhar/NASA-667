# -*- coding: utf-8 -*-
"""
Created on Thu Apr 06 13:10:28 2017

@author: aeber
"""
#simulates an observing campaign for various values of eta for a given telescope
import math
import random
import pickle
import Math
import rand_orbits as ro
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sp

#telescope stats
l, dl = .76e-6, .76e-7 #observed wavelength, and bandwidth
D = 8. #aperture diameter, m
scope = "HabEx" #needs to be JWST, WFIRST, HabEx, or LUVOIR
T = .1 #transmission
T_tot = 60.*60.*24.*365. #total mission lifetime, seconds
C = 1.0e-9 #contrast
IWA, OWA = 3., 20. #inner and outer working angles in l/D

#statistics information
z = 5. #target z score
eta = [.0001, .001, .01, .1, .3, .5, .7, .9] #true alue of eta
t_crit = 3.7 #gamma score value excluding at 95%
T_limit = T_tot*.17 #what is the maximum fraction of a mission we are willing to dedicate to one target
trials = 9999 #number of missions to run
t_steps_total = 300.

#prints some information about the simulation to help me keep track of multiple consoles
print scope
print T_limit/float(T_tot)
print t_crit
print t_steps_total

file_ = "hygdata_v3.csv" #file contain stellar info



#finds the fraction of the flux inside the bandwidth
#uses a blackbody approximation
def find_f(T):
    C = 2*math.pi*Math.h*dl*(Math.c**2)/float(Math.steph_boltz*(T**4)*(l**5))
    fac = math.exp(Math.h*Math.c/float(l*T*Math.k)) - 1.
    if fac == 0:
        fac = Math.h*Math.c/float(l*T*Math.k)
    return C/float(fac)

#loads the pickle containing information about stars
#this pickle file is made by find_exclude.py
def get_stars():
    stars = pickle.load(open("stars_" + scope + ".p", 'rb'))
    stars.sort()
    return stars

#simulates an observing campaign
#given star array (stars), eta value (eta_), discovered planet luminosity array (L),
#discovered planet semi major axis array (a), discovered planet period array (P),
#discovered planet z score array (z_scores), and arrays for the total planets created for each of those variables
#returns number of planets found and proportion found
#also edits the arrays passed to it
def runMission(stars, eta_, L, d, a, P, z_scores, L_tot, d_tot, a_tot, P_tot):
    found  = 0 #planets detected
    created = 0 #planets created
    T_obs = 0 #total time already observed
    for star_ in stars: #star loop
        if T_obs < T_tot: #checks if there is still time in mission
            if random.uniform(0,1) < eta_: #check to create planet
                created += 1
                L_tot.append(star_.L/float(Math.L_sun))
                d_tot.append(star_.d/float(Math.pc))
                a_tot.append(star_.planet.a/float(Math.au))
                P_tot.append(star_.planet.T/float(60.*60.*24.*365.))
                star_.add_planet()
                signal = 0 #signal total
                background = 0 #epextation for background
                t_step = star_.tau_exclude*T_tot*t_crit/float(t_steps_total)
                z_, looked = 0, 0
                while(looked < t_steps_total and z_ < z and t_step*looked < T_limit): #observing loop
                    light = star_.get_photons(D, C, l, IWA, OWA, t_step, T)
                    background += light[0]
                    signal += light[3] + light[2]
                    z_ = (signal - background)/float(math.sqrt(background))
                    star_.planet.progress(t_step)
                    looked += 1
                    if z_ > z: #if a planet is detected
                        found += 1
                        L.append(star_.L/float(Math.L_sun))
                        d.append(star_.d/float(Math.pc))
                        a.append(star_.planet.a/float(Math.au))
                        P.append(star_.planet.T/float(60.*60.*24.*365.))
                z_scores.append(z_)
                factor = 1.
                if z_ > 0: 
                    factor = min(1., (z/float(z_))**2) #corrects for over spent time
                #T_obs += looked*star_.tau*t_crit/float(t_steps_total)
                T_obs += factor*looked*star_.tau_exclude*T_tot*t_crit/float(t_steps_total)
            else: #if no planet is created in the first place
                T_obs += t_crit*star_.tau_exclude*T_tot
    try:
        return found, found/float(created)
    except ZeroDivisionError: #if no planets were created
        return 0, 0

#runs an MC sim for this eta value
#given star array (stars), eta value (eta_), discovered planet luminosity array (L),
#discovered planet semi major axis array (a), discovered planet period array (P),
#discovered planet z score array (z_scores), and arrays for the total planets created for each of those variables
#returns the average number found and average proportion found
def run_MC(stars, eta_, L, d, a, P, z_scores, L_tot, d_tot, a_tot, P_tot):
    found, proportion_found = 0, []
    for i in range(trials): #run trial missions
        results = runMission(stars, eta_, L, d, a, P, z_scores, L_tot, d_tot, a_tot, P_tot)
        found += results[0]
        proportion_found.append(results[1])
    return found/float(trials), np.mean(proportion_found)
    

#calculates a confidence interval for an estimator of eta given the simulation results
#returns an array with the estimator values, the upper limits on those values and the lower limits
def eta_est(expect):
    eta_hat, high, low = [], [], []
    for i in range(len(expect)):
        found = expect[i]        
        eta_ = eta[i]
        if found > 0:
            sig = math.sqrt(found)
            t_star = sp.t.ppf(.95, found)
            found_up = t_star*sig + found
            found_down = found - t_star*sig 
            high.append(min(eta_*found_up/float(found), 1.))
            low.append(max(0.0, eta_*found_down/float(found)))
            eta_hat.append(eta_)
        else:
            eta_hat.append(0)
            high.append(0)
            low.append(0)
    return eta_hat, high, low
    

#just divides each element in array A by the corresponding element in array B
#returns an array containing the result
def divide(A, B):
    C = []
    for i in range(len(A)):
        if B[i] != 0:
            C.append(A[i]/float(B[i]))
        else:
            C.append(0)
    return C

#returns an estimate of the transit probability for each value in an array of log semi major axis values
def getTP(logA):
    TP = []
    for i in range(len(logA)):
        a = Math.au*(10**logA[i])
        TP.append((Math.r_sun/float(a))**2)
    return TP


if __name__ == "__main__":
    stars = get_stars()
    planet_L, planet_d, planet_a, planet_P, find_expect = [], [], [], [], []
    L_tot, d_tot, a_tot, P_tot = [], [], [], []
    sensitivity, z_scores = [], []
    for eta_ in eta: #eta loop
        print "running sim for eta = ", eta_
        results = run_MC(stars, eta_, planet_L, planet_d, planet_a, planet_P, z_scores, L_tot, d_tot, a_tot, P_tot)
        find_expect.append(results[0])
        sensitivity.append(results[1])
        
    plt.figure(1)
    plt.title(scope + str(T_limit/float(T_tot)))
    plt.xlabel("Luminosity (log solar luminosity)")
    plt.ylabel("counts")
    A_L = plt.hist(np.log10(L_tot))
    B = plt.hist(np.log10(planet_L), bins = A_L[1])
    sensitivity_L = divide(B[0], A_L[0])
    plt.show()
    
    plt.figure(2)
    plt.title(scope + str(T_limit/float(T_tot)))
    plt.xlabel("Distance (pc)")
    A_d = plt.hist(d_tot)
    B = plt.hist(planet_d, bins = A_d[1])
    sensitivity_d = divide(B[0], A_d[0])
    plt.show()
    
    plt.figure(3)
    plt.title(scope + str(T_limit/float(T_tot)))
    plt.xlabel("log solar luminosity")
    plt.ylabel("sensitivity")
    plt.plot(A_L[1][1:], sensitivity_L)
    print "max found: ", max(find_expect)
    plt.show()
    
    plt.figure(4)
    plt.title(scope + str(T_limit/float(T_tot)))
    plt.xlabel("distance (pc)")
    plt.ylabel("sensitivity")
    plt.plot(A_d[1][1:], sensitivity_d)
    print "sensitivity: ", np.median(sensitivity)
    plt.show()
    
    plt.figure(5)
    plt.title(scope + str(T_limit/float(T_tot)))
    plt.xlabel("z scores")
    plt.hist(z_scores, bins = 100)
    plt.show()
    
    plt.figure(6)
    plt.title(scope + str(T_limit/float(T_tot)))
    plt.xlabel("log solar luminosities")
    plt.ylabel("distances (pc)")
    plt.scatter(np.log10(planet_L), planet_d)
    plt.show()
    
    eta_hat, high, low = eta_est(find_expect)
    plt.figure(7)
    plt.xlabel("true value of eta")
    plt.ylabel("measured eta")
    plt.plot(eta, eta_hat)
    plt.plot(eta, high)
    plt.plot(eta, low)
    plt.show()
    
    plt.figure(8)
    plt.xlabel("semi major axis (au)")
    A_a = plt.hist(np.log10(a_tot), bins = 100)
    B = plt.hist(np.log10(planet_a), bins = A_a[1])
    sensitivity_a = divide(B[0], A_a[0])
    plt.show()
    
    plt.figure(9)
    plt.xlabel("log period (yrs)")
    A_P = plt.hist(np.log10(P_tot), bins = 100)
    B = plt.hist(np.log10(planet_P), bins = A_P[1])
    sensitivity_P = divide(B[0], A_P[0])
    plt.show()
    
    transit_P = getTP(A_a[1][1:])
    
    plt.figure(10)
    plt.xlabel("log semi major axis (au)")
    plt.ylabel("sensitivity")
    plt.plot(A_a[1][1:], sensitivity_a)
    plt.plot(A_a[1][1:], transit_P)
    plt.show()
    
    plt.figure(11)
    plt.xlabel("log period (yrs)")
    plt.ylabel("sensitivity")
    plt.plot(A_P[1][1:], sensitivity_P)
    plt.show()
