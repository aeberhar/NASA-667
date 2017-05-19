# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 10:06:58 2016

@author: aeber
"""
import candidate1_1 as candidate
import manipulateAll as mA
import matplotlib.pyplot as plt
direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys
sys.path.insert(0, direct)
import Math

def modelStar(T_eff, r_star, L_file, D):
    o = candidate.ToI()
    o.T_eff = T_eff
    o.r_star = r_star
    o.L_file = L_file
    o.D = D
    return o


if __name__ == "__main__":
    l = range(20, 4000, 10)
    dl = 1e-8
    plt.figure(1)
    Barnard = modelStar([3134,102], [.196*Math.r_solar, .008*Math.r_solar], "kp00_3500.ascii", [1.8328*Math.pc, .0007*Math.pc])
    L_Barnard = mA.findStellarSpectrum(Barnard, l, dl, True)
    print "Barnard's star:", mA.getHZ(Barnard)/Math.au, float(mA.getHZ(Barnard)/Barnard.D[0])/(.5e-6 / 6.), (Math.r_earth/mA.getHZ(Barnard))**2
    plt.plot(l, L_Barnard, color = 'b')
    SiriusA = modelStar([9940,0], [1.711*Math.r_solar, .0*Math.r_solar], "kp00_9750.ascii", [2.64*Math.pc, .01*Math.pc])
    L_SiriusA = mA.findStellarSpectrum(SiriusA, l, dl, True)
    print "Sirius:", mA.getHZ(SiriusA)/Math.au, float(mA.getHZ(SiriusA)/SiriusA.D[0])/(.5e-6 / 6.), (Math.r_earth/mA.getHZ(SiriusA))**2
    plt.plot(l, L_SiriusA, color = 'g')
    SiriusB = modelStar([25200,0], [.0084*Math.r_solar, .0*Math.r_solar], "kp00_26000.ascii", [2.64*Math.pc, .01*Math.pc])
    L_SiriusB = mA.findStellarSpectrum(SiriusB, l, dl, True)
    plt.plot(l, L_SiriusB, color = 'm')
    EpsEri = modelStar([5084,6], [.735*Math.r_solar, .005*Math.r_solar], "kp00_5250.ascii", [3.212*Math.pc, .001*Math.pc])
    L_EpsEri = mA.findStellarSpectrum(EpsEri, l, dl, True)
    print "EpsEri:", mA.getHZ(EpsEri)/Math.au, float(mA.getHZ(EpsEri)/EpsEri.D[0])/(.5e-6 / 6.), (Math.r_earth/mA.getHZ(EpsEri))**2
    plt.plot(l, L_EpsEri, color = 'r')
    ProcyonA = modelStar([6530,50], [2.048*Math.r_solar, .025*Math.r_solar], "kp00_6500.ascii", [3.51*Math.pc, .02*Math.pc])
    L_ProcyonA = mA.findStellarSpectrum(ProcyonA, l, dl, True)
    print "Procyon:", mA.getHZ(ProcyonA)/Math.au, float(mA.getHZ(ProcyonA)/ProcyonA.D[0])/(.5e-6 / 6.), (Math.r_earth/mA.getHZ(ProcyonA))**2
    plt.plot(l, L_ProcyonA, color = 'c')
    ProcyonB = modelStar([7740,50], [.01234*Math.r_solar, .00032*Math.r_solar], "kp00_7750.ascii", [3.51*Math.pc, .02*Math.pc])
    L_ProcyonB = mA.findStellarSpectrum(ProcyonB, l, dl, True)
    plt.plot(l, L_ProcyonB, color = 'y')
    TauCeti = modelStar([5344,50], [.793*Math.r_solar, .004*Math.r_solar], "kp00_5250.ascii", [3.650*Math.pc, .002*Math.pc])
    L_TauCeti = mA.findStellarSpectrum(TauCeti, l, dl, True)
    print "TauCeti:", mA.getHZ(TauCeti)/Math.au, float(mA.getHZ(TauCeti)/TauCeti.D[0])/(.5e-6 / 6.), (Math.r_earth/mA.getHZ(TauCeti))**2
    plt.plot(l, L_TauCeti, color = 'k')
    Groombridge = modelStar([3970,0], [.605*Math.r_solar, .002*Math.r_solar], "kp00_4000.ascii", [3.587*Math.pc, .010*Math.pc])
    L_Groombridge = mA.findStellarSpectrum(Groombridge, l, dl, True)
    print "Groombridge:", mA.getHZ(Groombridge)/Math.au, float(mA.getHZ(Groombridge)/Groombridge.D[0])/(.5e-6 / 6.), (Math.r_earth/mA.getHZ(Groombridge))**2
    plt.plot(l, L_Groombridge, color = 'k')    
    plt.show()
    plt.figure(2)
    Math.scaleMatrix(L_Barnard, 1./(4*Math.pi*Barnard.D[0]**2))
    plt.plot(l, L_Barnard)
    Math.scaleMatrix(L_SiriusA, 1./(4*Math.pi*SiriusA.D[0]**2))
    plt.plot(l, L_SiriusA)
    Math.scaleMatrix(L_SiriusB, 1./(4*Math.pi*SiriusB.D[0]**2))
    plt.plot(l, L_SiriusB)
    Math.scaleMatrix(L_EpsEri, 1./(4*Math.pi*EpsEri.D[0]**2))
    plt.plot(l, L_EpsEri)
    Math.scaleMatrix(L_ProcyonA, 1./(4*Math.pi*ProcyonA.D[0]**2))
    plt.plot(l, L_ProcyonA)
    Math.scaleMatrix(L_ProcyonB, 1./(4*Math.pi*ProcyonB.D[0]**2))
    plt.plot(l, L_ProcyonB)
    Math.scaleMatrix(L_TauCeti, 1./(4*Math.pi*TauCeti.D[0]**2))
    plt.plot(l, L_TauCeti)
    Math.scaleMatrix(L_Groombridge, 1./(4*Math.pi*Groombridge.D[0]**2))
    plt.plot(l, L_Groombridge)    
    plt.show()