# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:43:19 2020
Age-specific SEIRD model.
@author: Jason
"""
import numpy as np

def model(EpidemicParameters, regionParameters):
    R0 = EpidemicParameters['R0']
    numAge = 17
    infectiousRateList = EpidemicParameters['infectiousRateList']
    confirmRateList = EpidemicParameters['confirmRateList']
    asymRecoverRateList = EpidemicParameters['asymRecoverRateList']
    mildRecoverRateList = EpidemicParameters['mildRecoverRateList']
    severeRecoveRateList = EpidemicParameters['severeRecoveRateList']
    mildToSevereRateList = EpidemicParameters['mildToSevereRateList']
    mortalityRateList = EpidemicParameters['mortalityRateList']
    initialState = regionParameters['initialState']
    population = regionParameters['population']
    populationRate = regionParameters['populationRate']
    CMList = regionParameters['contactMatrixList']
    normalizeCM = regionParameters['normalizeContactMatrix']
    duration = regionParameters['duration']
    
    def getBeta(R0, CM, recoveryRate, calTransProb=True):
        if calTransProb:
            M = CM
            eigenvalue, featurevector = np.linalg.eig(np.matrix(M))
            beta = R0 * recoveryRate / max(np.real(eigenvalue))
        else:
            beta = 0.025
        return beta
    
    def exposedSum(age, prevSimuResult, beta, CM, asymptomatic):
        sumprevSimuValue = 0
        if asymptomatic:
            for j in range(numAge):
                sumprevSimuValue += beta * CM[age,j] * prevSimuResult[j*7+2] * prevSimuResult[age*7+0] / population[j]
        else:
            for j in range(numAge):
                sumprevSimuValue += beta * CM[age,j] * prevSimuResult[j*7+3] * prevSimuResult[age*7+0] / population[j]
        return sumprevSimuValue

    betaAsym = getBeta(R0, sum(normalizeCM), asymRecoverRateList[0])
    #beta_s = get_beta(R0, sum(region.normalize_cm), mildRecoverRateList[0], populationRate)
    betaMild = betaAsym

    def seirmicr(prevSimuResult, dt, t):
        curSimuResult = np.zeros((numAge*7))
        dailyCase = np.zeros((numAge*7))      
        CM = np.zeros((numAge,numAge))
        try:
            CM = CMList[int(t)]
        except:
            print('Cannot achieve contact matrix!')
            print(int(t))
            
        for age in range(numAge):
            asymExposed = exposedSum(age, prevSimuResult, betaAsym, CM, asymptomatic=True) * dt
            mildExposed = exposedSum(age, prevSimuResult, betaMild, CM, asymptomatic=False) * dt
            exposed = asymExposed + mildExposed
            infection = infectiousRateList[age]*prevSimuResult[age*7+1] *dt
            asymInfection = (1-confirmRateList[age]) * infection
            mildInfection = confirmRateList[age] * infection
            asymRecover = asymRecoverRateList[age] * prevSimuResult[age*7+2] * dt
            mildRecover = mildRecoverRateList[age] * (1-mildToSevereRateList[age]) * prevSimuResult[age*7+3] * dt
            mildToSevereInfection = mildToSevereRateList[age] * prevSimuResult[age*7+3] * dt
            severeRecover = severeRecoveRateList[age] * (1-mortalityRateList[age]) * prevSimuResult[age*7+4] * dt
            mildDeath = mortalityRateList[age] * prevSimuResult[age*7+3] *dt
            severeDeath = mortalityRateList[age] * prevSimuResult[age*7+4] * dt
            
            dailyCase[age*7+0] = -exposed
            dailyCase[age*7+1] = exposed
            dailyCase[age*7+2] = asymInfection
            dailyCase[age*7+3] = mildInfection
            dailyCase[age*7+4] = mildToSevereInfection
            dailyCase[age*7+5] = asymRecover + mildRecover + severeRecover
            dailyCase[age*7+6] = mildDeath + severeDeath

            curSimuResult[age*7+0] = max(prevSimuResult[age*7+0] - exposed, 0)
            curSimuResult[age*7+1] = prevSimuResult[age*7+1] + exposed - infection
            curSimuResult[age*7+2] = prevSimuResult[age*7+2] + asymInfection - asymRecover
            curSimuResult[age*7+3] = prevSimuResult[age*7+3] + mildInfection - mildRecover - mildToSevereInfection - mildDeath
            curSimuResult[age*7+4] = prevSimuResult[age*7+4] + mildToSevereInfection - severeRecover - severeDeath
            curSimuResult[age*7+5] = prevSimuResult[age*7+5] + asymRecover + mildRecover + severeRecover
            curSimuResult[age*7+6] = prevSimuResult[age*7+6] + mildDeath + severeDeath
        return curSimuResult, dailyCase
    
    def splitResult(result, timeIdx, u):
        for age in range(numAge):
            result[0][timeIdx, age] = u[age*7+0]
            result[1][timeIdx, age] = u[age*7+1]
            result[2][timeIdx, age] = u[age*7+2]
            result[3][timeIdx, age] = u[age*7+3]
            result[4][timeIdx, age] = u[age*7+4]
            result[5][timeIdx, age] = u[age*7+5]
            result[6][timeIdx, age] = u[age*7+6]
        return result

    u = initialState
    dt = 1 
    tl = int((duration-1)*(1/dt))+1
    t = np.linspace(0, duration-1, tl)

    S = np.zeros((tl, numAge))
    E = np.zeros((tl, numAge))
    Iasym = np.zeros((tl, numAge))
    Im = np.zeros((tl, numAge))
    Is = np.zeros((tl, numAge))
    R = np.zeros((tl, numAge))
    D = np.zeros((tl, numAge))
    cumCase = [S, E, Iasym, Im, Is, R, D]
    cumCase = splitResult(cumCase, 0, u)
         
    d_s = np.zeros((tl, numAge))
    d_e = np.zeros((tl, numAge))
    d_iasym = np.zeros((tl, numAge))
    d_im = np.zeros((tl, numAge))
    d_is = np.zeros((tl, numAge))
    d_r = np.zeros((tl, numAge))
    d_d = np.zeros((tl, numAge)) 
    dailyCase = [d_s, d_e, d_iasym, d_im, d_is, d_r, d_d]

    for j in range(1, tl):
        u, daily = seirmicr(u, dt, t[j])
        cumCase = splitResult(cumCase, j, u)
        dailyCase = splitResult(dailyCase, j, daily)

    return cumCase, dailyCase
