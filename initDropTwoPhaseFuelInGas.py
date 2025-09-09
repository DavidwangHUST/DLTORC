#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 14:39:59 2024

@author: wangweiye

For binary component,we allow fuel species exist in gas phase
"""

import numpy as np
from numpy import power
import cantera as ct
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline

font = {'family':'times',
        'color':'darkred',
        'weight':'normal',
        'size':16}
#%%

def calVaporPressure(T:float) :
    #    return 1e5*np.power(10,4.6543 - (1435.264/(T - 64.84))) # returns Pa, where T is in K
    #P = 1e3*np.exp(16.7 - (4060.0/(T - 37.0))) # returns Pa, where T is in K,FOR WATER
    #P = 1e5*np.power(10,4.02832-(1268.636/(T-56.199))) ; #FOR N-HEPTANE
    propaneVapPres = 1.0e5*pow(10,4.53678-(1149.36/(T+24.906))) # propane vapor equation 
    heptaneVapPres = 1.0e5*pow(10,4.02832-(1268.636/(T-56.199))) # heptane vapor equation 
    octaneVapPres = 1.0e5*pow(10,3.93679-(1257.84/(T-52.415)))
    # propaneVapPres = 0.00  #NO FUEL @ t0
    # heptaneVapPres = 0.00 

    return propaneVapPres,heptaneVapPres

# return 2D array which contains fuel and "air" mole fractions at left and right boundary
def genGasFuelAmbientMoleArr(T,P,gas,gasComp,fuelComp):
    gas.TPX = T,P, fuelComp
    propaneMoleFrac,heptaneMoleFrac = calVaporPressure(T)
    fuelMoleArray = np.array([gas[fuel[0]].X[0] , gas[fuel[1]].X[0]])

    propaneMoleFrac = fuelMoleArray[0]/fuelMoleArray.sum() * propaneMoleFrac/ P
    heptaneMoleFrac = fuelMoleArray[1]/fuelMoleArray.sum() * heptaneMoleFrac/ P

    weightFactor = 1.00 - (propaneMoleFrac + heptaneMoleFrac)
    gas.TPX = T,P, gasComp
    diluentMoleSum = gas['O2'].X[0] + gas['N2'].X[0]
    oxygenMoleFracRight =  gas['O2'].X[0] / diluentMoleSum
    nitrogenMoleFracRight = 1.00 - oxygenMoleFracRight
    oxygenMoleFracLeft = oxygenMoleFracRight * weightFactor
    nitrogenMoleFracLeft = nitrogenMoleFracRight * weightFactor

    speciesMoleFracArr = np.array([[propaneMoleFrac, 0],
                                   [heptaneMoleFrac, 0],
                                   [oxygenMoleFracLeft, oxygenMoleFracRight],
                                   [nitrogenMoleFracLeft,nitrogenMoleFracRight]])
    return speciesMoleFracArr

# return 2D full array with 
def genGasMassFracArr(Tdrop,Tgas,P,gas,gasComp,dropComp,Rd,L,shift,Wmix,nPts):
    dX = L/(nPts-1)
    r_ = np.zeros(nPts)
    r_[0] = Rd
    #DEBUG
    #print("r_[0] = %15.6e\n"%(r_[0]))
    for ii in range(1, nPts-1):
        r_[ii] = r_[ii-1] + dX
    r_[nPts-1] = Rd+L

    speciesMoleFracArr = genGasFuelAmbientMoleArr(Tdrop, P, gas, gasComp,dropComp)
    gasMoleArray = np.zeros([gas.n_species,nPts])
    gasMassArray = np.zeros([gas.n_species,nPts])
    gas.TPX = Tgas,P, dropComp +','+gasComp
    targetSpeciesIndexArray = np.array([gas.species_index(fuel[0]),
                                        gas.species_index(fuel[1]),
                                        gas.species_index('O2'),
                                        gas.species_index('N2')])

    for rowII in range(speciesMoleFracArr.shape[0]):
        x = [Rd+shift, Rd+shift+Wmix]
        y = [speciesMoleFracArr[rowII,:][0], speciesMoleFracArr[rowII,:][1]]
        m = [0,0]
        spline = CubicHermiteSpline(x, y, m)

        tempY_ = np.zeros(nPts)
        for i in range(nPts):
            if r_[i] < Rd + shift :
                tempY_[i] = y[0] ;
            if r_[i] >= (Rd + shift) and (r_[i] <= Rd + shift + Wmix) :
                tempY_[i] = spline(r_[i]) ;
            if r_[i] > (Rd + shift +Wmix) and r_[i] <= (Rd+L) :
                tempY_[i] = y[1] ;
        gasMoleArray[targetSpeciesIndexArray[rowII],:] = tempY_[:]

    for colummII in range(nPts):
        gas.TPX = Tgas,P,gasMoleArray[:,colummII]
        gasMassArray[:,colummII] = gas.Y

    # #DEBUG
    # fig,ax = plt.subplots(nrows=2,ncols=1,dpi=200)
    # ax[0].plot(np.array(r_/Rd),gasMassArray[targetSpeciesIndexArray[1],:],lw=2.0)
    # ax[0].set_ylabel("heptane mass fraction")
    # ax[0].set_xlim([0.0,8.0])
    # ax[1].plot(r_/Rd,gasMassArray[targetSpeciesIndexArray[-1],:],lw=2.0)
    # ax[1].set_xlim([0.0,8.0])
    # ax[1].set_ylabel("nitrogen mass fraction")
    # plt.show()

    return gasMassArray


# return 1D numpy array for liquid phase mass fraction
def genLiquidMassFracArr(gas,dropTemp:float,P:float,dropSpec:str) :
    gas.TPX = dropTemp,P,dropSpec
    massArr_ = gas.Y
    # print(gas.Y)
    return massArr_

# return 1D numpy array for gas phase mass fraction
# def genGasMassFracArr(gas,gasTemp:float,P:float,gasSpec:str):
#     gas.TPX = gasTemp,P,gasSpec
#     massArr_ = gas.Y
#     return massArr_

# return 1D numpy array for gas phase spatial coordinate and temperature
def genGasTempAndRadiusArr(Rd,L,shift,Wmix,nPts,dropTemp,gasTemp):
    dX = L/(nPts-1)
    r_ = np.zeros(nPts)
    r_[0] = Rd
    #DEBUG
    #print("r_[0] = %15.6e\n"%(r_[0]))
    for ii in range(1, nPts-1):
        r_[ii] = r_[ii-1] + dX
    r_[nPts-1] = Rd+L

    x = [Rd+shift, Rd+shift+Wmix] ;
    y = [dropTemp, gasTemp] ;
    m = [0,0] ;
    spline = CubicHermiteSpline(x, y, m)

    T_ = np.zeros(nPts)
    for i in range(nPts):
        if r_[i] < Rd + shift :
            T_[i] = dropTemp ;
        if r_[i] >= (Rd + shift) and (r_[i] <= Rd + shift + Wmix) :
            T_[i] = spline(r_[i]) ;
        if r_[i] > (Rd + shift +Wmix) and r_[i] <= (Rd+L) :
            T_[i] = gasTemp ;

    # #DEBUG
    # fig,ax = plt.subplots(dpi=200)
    # ax.plot(r_/Rd,T_,linewidth=2.0)
    # ax.set_xlim([0.0,8.0])
    # ax.set_xlabel("scaled radial coordinate")
    # ax.set_ylabel("temperature [K]")
    # plt.show()

    return r_,T_

# return 2D numpy array with size of (nvar*nPts)
def writeGlobalArr(Rd,L,shift,Wmix,dropTemp,gasTemp,P,dropSpec,gasSpec,gas,lnPts,gnPts,gasFileName,liquidFileName):
    liquidMassFracArr_ = genLiquidMassFracArr(gas, dropTemp, P, dropSpec)
    # gasMassFracArr_ = genGasMassFracArr(gas, gasTemp, P, gasSpec,dropSpec)
    gasMassFracArr_ = genGasMassFracArr(dropTemp,gasTemp, P, gas, gasSpec, dropSpec, Rd, L, shift, Wmix, gnPts)

    # assign spatial coordinate value
    rG_,TG_ = genGasTempAndRadiusArr(Rd, L, shift, Wmix, gnPts, dropTemp, gasTemp)
    rL_ = np.zeros(lnPts)
    dXL = Rd/(lnPts-1)
    rL_[0] = 0.0
    for ii in range(1, lnPts-1):
        rL_[ii] = rL_[ii-1] + dXL
    rL_[lnPts-1] = Rd
    #DEBUG
    # print("last element of rL_ is :%15.6e"%(rL_[lnPts-1]))
    # print("first element of rG_ is :%15.6e"%(rG_[0]))
    r_ = np.concatenate((rL_,rG_))
    # assign temperature values
    TL_ = np.ones(lnPts)*dropTemp
    T_ = np.concatenate((TL_,TG_))
    # assign pressure values
    P_ = np.ones(gnPts+lnPts) * P
    mdot_ = np.ones(gnPts+lnPts) * 0.0

    # write global arr to output file (both liquid phase and gas phase)
    out1=open(liquidFileName,"w")
    for i in range(lnPts):
        out1.write("%15.6e\t%15.6e\t"%(r_[i],T_[i]))
        for j in range(gas.n_species):
            out1.write("%15.6e\t"%(liquidMassFracArr_[j]))
        # out.write("%15.6e\t"%(P_[i]))
        out1.write("%15.6e\n"%(P_[i]))
    out1.close()
    # out.write("%15.6e\n"%(mdot_[i]))


    out2=open(gasFileName,'w')
    for i in range(lnPts,lnPts+gnPts):
        out2.write("%15.6e\t%15.6e\t"%(r_[i],T_[i]))
        for j in range(gas.n_species):
            out2.write("%15.6e\t"%(gasMassFracArr_[j,i-lnPts]))
        out2.write("%15.6e\n"%(P_[i]))
        # out.write("%15.6e\n"%(mdot_[i]))
    out2.close()


# =============================================================================
# Following block defines the user defined parameters
# =============================================================================
if __name__ == "__main__" :
    gasfileName = "initialCondition.dat"   #This script output the initial conditions
    liquidfileName="initialLiquidCondition.dat"
    #mech='sdMechMergedVer2.cti'
    mech='./NH3_NC7.xml'      #Mechanism name
    #mech='chem171.cti'
    gas = ct.Solution(mech)
    # gas1 = ct.Solution(mech)
    Rd = 50.0e-6                #droplet radius, Unit:[m]
    domainLength = Rd*50.0      #domain length, Unit:[m]
    shift = Rd * 0.0e0              #shift of mixing region from droplet surface
    Wmix = Rd * 0.5               #thickness of mixing layer for heptane mass fraction & Temperature distribution, Unit:[m]
    gNpts = 5000                #Number of grid points in gas phase
    lNpts = 500                 #Number of grid points in liquid phase
    Tamb = 850.0                #Temperature of ambient air, Unit : [K]
    Pamb = 20.0*ct.one_atm          #Pressure of ambient air, Unit :[Pa]
    Td = 450.00                 #Temperature of droplet surface, Unit: [K]
    #RH = 0.0                    #relative humidity
    #dropName = "C3H8:0.20,NC7H16:0.80"          #Name of droplet species
    fuel = ["C3H6", "NC7H16"]
    fuelMole = [0.00,1.00]
    dropName = ",".join(f"{f}:{m}" for f, m in zip(fuel, fuelMole))
    # print(result)  # Output: "C3H8:0.10,NC7H16:0.90"
    #dropName = "NC7H16:1.00"
    #gasComposition = "O2:0.21,N2:0.79"
    gasComposition = "NO2:0.1,O2:20.586,N2:77.443"

    writeGlobalArr(Rd, domainLength, shift, Wmix, Td, Tamb, Pamb, dropName, gasComposition, gas, lNpts, gNpts, gasfileName,liquidfileName)
    # writeGlobalArr(Rd, L, shift, Wmix, dropTemp, gasTemp, P, dropSpec, gasSpec, gas, lnPts, gnPts, fileName)
    #writeInitCond(Rd, domainLength, shift, Wmix, nGrid, Tamb, Pamb, Td, RH, gas, dropName, oxidizerName, bathName, fileName)
    # genTotalMassfracArr(Rd, domainLength, shift, Wmix, nGrid, Tamb, Pamb, Td, RH, gas, dropName, oxidizerName, bathName)
    # genTemArr(Rd, domainLength, shift, Wmix, nGrid, Tdrop, Tamb)
    # genMolefracArr(Rd, domainLength, shift, Wmix, nGrid, Tamb, Pamb, Td, RH)
