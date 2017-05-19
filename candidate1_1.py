# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:53:06 2016

@author: aeber
"""


direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys
sys.path.insert(0, direct)
import Math

class ToI(object):
    
    def __init__(self, earth = False):
        self.n = 0 #candidate number
        self.p = [] #p values for arbitrary cuts
        self.d = 10 #width of apeture in meters
        self.R = 100 #spectral resolution l/dl
        self.S = 1e-5
        self.T = 100 #observation time in s
        self.N_str = [["N2", 2.0e22], ["O3", 1.5e20], ["O2", 4.0e21]] #gas content matrix for stratosphere
        self.t_str = 30.0e3 #thickness of statosphere in meters
        self.N_tr= [["N2", 2.0e25], ["O2", 5.0e24], ["H2O", 1.0e22]]
        self.t_tr = 20.0e3
        if earth:
            self.r_p = [Math.r_earth,0] #planet radius
            self.L_file = "kp00_5750.ascii" #file for kurusc data
            self.a = [Math.au, 0] #orbital radius
            self.m = [Math.m_earth,0] #planet mass
            self.T_eff = [5770,0] #effective temperature of star
            self.r_star = [Math.r_solar,0] #radius of star
            self.D = [10.0*Math.pc, 0] #distance to system
        else:
            self.r_p = [0,0] #planet radius
            self.a = [0, 0] #orbital radius
            self.m = [0,0] #planet mass
            self.T_eff = [0,0] #effective temperature of star
            self.L_file = False
            self.r_star = [0,0] #radius of star
            self.D = [0, 0] #distance to system
    
    def __str__(self):
        print "planet mass: ", self.m
    
    #copies the attributes of given object o to this object
    def copy(self,o):
        self.n = o.n
        self.r_p = o.r_p
        self.a = o.a
        self.m = o.m
        self.T_eff = o.T_eff
        self.r_star = o.r_star
        self.D = o.D
        self.p = o.p
        self.N_str = o.N_str
        self.t_str = o.t_str
        self.N_tr = o.N_tr
        self.t_tr = o.t_tr
        self.d = o.d


#class for candidate object
class candidate(ToI):
    
    def __init__(self, n):
        self.n = n
        self.L = [0,0] #luminosity of star
        self.mag = [0,0] #magnitude of star
    
    #tries to find luminosity
    def getL(self):
        if self.T_eff[0] and self.r_star[0]:
            T4, sT4 = Math.square(self.T_eff[0], self.T_eff[1], 4)
            r2, sr2 = Math.square(self.r_star[0], self.r_star[1])
            T4r2, sT4r2 = Math.multiply(T4, sT4, r2, sr2)
            L, sL = Math.scale(T4r2, sT4r2, 4*Math.math.pi*Math.steph_boltz)
            if not(self.L[0]) or (self.L[1] > sL):            
                self.L[0], self.L[1] = L, sL
        if self.D[0] and self.mag[0]:
            alpha, salpha = (self.mag[0] + 27)/-2.5, abs(self.mag[1]/2.5)
            D2, sD2 = Math.square(self.D[0], self.D[1])
            _10a, s10a = Math.pow(10, 0, alpha, salpha)
            D210a, sD210a = Math.multiply(D2, sD2, _10a, s10a)
            L, sL = Math.scale(D210a, sD210a, Math.L_solar/(Math.au**2))
            if not(self.L[0]) or (self.L[1] > sL):            
                self.L[0], self.L[1] = L, sL
           
    
    #tries to find distance to system   
    def getD(self):
        if self.L[0] and self.mag[0]:
            alpha, salpha = (self.mag[0] + 27)/-2.5, abs(self.mag[1]/2.5)
            _10a, s10a = Math.pow(10, 0, alpha, salpha)
            L_10a, sL_10a = Math.divide(self.L[0], self.L[1], _10a, s10a)
            _radical, s_radial = Math.scale(L_10a, sL_10a, (Math.au**2)/Math.L_solar)
            D, sD = Math.sqrt(_radical, s_radial)
            if not(self.D[0]):
                self.D[0], self.D[1] = D, sD
    
    
    #tries to find the orbital semimajor axis    
    def geta(self):
        if not(self.a[0]) and self.D[0]:
            self.a[0], self.a[1] = Math.scale(self.D[0], self.D[1], 2.0*(.5e-6)/6.0)
    
    #fills in remaining necessary planet info
    def remainder(self):
        if not(self.r_p[0]):
            self.r_p[0] = [Math.r_earth,0]
        if not(self.m[0]):
            self.m[0] = [Math.m_earth,0]
        if not(self.a[0]):
            self.a[0] = [Math.au,0]
        if not(self.T_eff[0]):
            self.T_eff[0] = [5750,0]
        if not(self.r_star[0]):
            self.r_star[0] = [Math.r_solar,0]
        if not(self.D[0]):
            self.D[0] = [Math.pc*2.0,0]
            