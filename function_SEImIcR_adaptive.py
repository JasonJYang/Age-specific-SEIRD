# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:43:19 2020
Age-specific SEIRD model for adaptive policy.
@author: Jason
"""
import numpy as np

def model(EpidemicParameters, regionParameters, adaptiveControl, Region):
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
    normalizeCM = regionParameters['normalizeContactMatrix']
    duration = regionParameters['duration']
    
    def getBeta(R0, CM, recoveryRate, calTransProb=True):
        if calTransProb:
            # get the eigenvalue and featurevector for the matrix M
            M = CM
            eigenvalue, featurevector = np.linalg.eig(np.matrix(M))
            beta = R0 * recoveryRate / max(np.real(eigenvalue))
        else:
            beta = 0.025
        return beta
    
    def generateCM(dailyCaseAll, t, controlTime):
        if 'severe' in adaptiveControl:
            infections = sum(dailyCaseAll[4][int(t)-1])
            threshold = adaptiveControl['severe']
        elif 'mild' in adaptiveControl:
            infections = sum(dailyCaseAll[3][int(t)-1])
            threshold = adaptiveControl['mild']
        elif 'death' in adaptiveControl:
            infections = sum(dailyCaseAll[6][int(t)-1])
            threshold = adaptiveControl['death']
        if infections >= threshold:
            weightDict = {'home': np.diag(list(np.repeat(1, numAge))),
                      'work': 0,
                      'school': 0,
                      'other': np.diag(list(np.repeat(0.1, numAge)))}
            controlTime.append(int(t))
        else:
            weightDict = {'home': 1,
                      'work': 1,
                      'school': 1,
                      'other': 1}
        return Region.getWeightedCM(weightDict)
    
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
    betaMild = betaAsym
    
    def seirmicr(prevSimuResult, dt, t, dailyCaseAll, controlTime):
        curSimuResult = np.zeros((numAge*7))
        dailyCase = np.zeros((numAge*7))      
        CM = generateCM(dailyCaseAll, t, controlTime)
            
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
    duration = duration - 1
    tl = int(duration*(1/dt))+1
    t = np.linspace(0, duration, tl)

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
    dailyCaseAll = [d_s, d_e, d_iasym, d_im, d_is, d_r, d_d]

    controlTime = []
    
    for j in range(1, tl):
        u, daily = seirmicr(u, dt, t[j], dailyCaseAll, controlTime)
        cumCase = splitResult(cumCase, j, u)
        dailyCaseAll = splitResult(dailyCaseAll, j, daily)

    return cumCase, dailyCaseAll, controlTime

def seir(covid, region, control):
    R0 = covid.R0
    n_age = covid.n_age
    #transmission_rate = covid.transmission_rate
    infectious_rate = covid.infectious_rate
    confirmed_rate = covid.confirmed_rate
    as_recovery_rate = covid.as_recovery_rate
    m_recovery_rate = covid.m_recovery_rate
    s_recovery_rate = covid.s_recovery_rate
    #severe_rate = covid.severe_rate
    ms_rate = covid.ms_rate
    motality_rate = covid.motality_rate
    initial_state = region.initial_state
    population = region.population
    population_rate = region.population_rate
    duration = region.duration
    P = population
    
    def get_beta(R0, contact_matrix, recovery_rate, 
                 population_rate, calculate_transmission_probability=True):
        if calculate_transmission_probability:
            M = contact_matrix
            for i in range(n_age):
                for j in range(n_age):
                    M[i,j] = contact_matrix[i,j] * population_rate[i] / population_rate[j]
            # get the eigenvalue and featurevector for the matrix M
            eigenvalue, featurevector = np.linalg.eig(np.matrix(M))
            beta = R0 * recovery_rate / max(np.real(eigenvalue))
        else:
            beta = 0.025
        return beta
    
    def si_sum(i, p, V, beta, cm, asi):
        sum_value = 0
        if asi:
            for j in range(n_age):
                sum_value += beta * cm[i,j] * V[j*7+2] * V[i*7+0] / p[j]
        else:
            for j in range(n_age):
                sum_value += beta * cm[i,j] * V[j*7+3] * V[i*7+0] / p[j]
        """
        if sum_value >= V[i*7+0]:
            sum_value = V[i*7+0]
        """
        return sum_value
    
    def generate_cm(daily_all, t, control_time):
        if 'severe' in control:
            infections = sum(daily_all[4][int(t)-1])
            threshold = control['severe']
        elif 'mild' in control:
            infections = sum(daily_all[3][int(t)-1])
            threshold = control['mild']
        elif 'death' in control:
            infections = sum(daily_all[6][int(t)-1])
            threshold = control['death']
        if infections >= threshold:
            weight = {'home': np.diag(list(np.repeat(1, n_age))),
                      'work': 0,
                      'school': 0,
                      'other': np.diag(list(np.repeat(0.1, n_age)))}
            control_time.append(int(t))
            return region.get_cm(weight)
        else:
            weight = {'home': 1,
                      'work': 1,
                      'school': 1,
                      'other': 1}
            return region.get_cm(weight)

    beta_as = get_beta(R0, sum(region.normalize_cm), as_recovery_rate[0], population_rate)
    #beta_s = get_beta(R0, sum(region.normalize_cm), m_recovery_rate[0], population_rate)
    beta_s = beta_as
    
    def seirmicr(V, dt, t, daily_all, control_time):
        Y = np.zeros((n_age*7))
        daily_case = np.zeros((n_age*7))
        
        cm = generate_cm(daily_all, t, control_time)
            
        for age in range(n_age):
            as_exposed = si_sum(age, P, V, beta_as, cm, asi=True) * dt
            s_exposed = si_sum(age, P, V, beta_s, cm, asi=False) * dt
            exposed = as_exposed + s_exposed
            infection = infectious_rate[age]*V[age*7+1] *dt
            as_infection = (1-confirmed_rate[age]) * infection
            m_infection = confirmed_rate[age] * infection
            as_recovery = as_recovery_rate[age] * V[age*7+2] * dt
            m_recovery = m_recovery_rate[age] * (1-ms_rate[age]) * V[age*7+3] * dt
            mtos_infection = ms_rate[age] * V[age*7+3] * dt
            #s_infection = severe_rate[age] * m_infection
            s_recovery = s_recovery_rate[age] * (1-motality_rate[age]) * V[age*7+4] * dt
            #s_death = motality_rate[age] * V[age*7+4] * dt
            m_death = motality_rate[age] * V[age*7+3] *dt
            s_death = motality_rate[age] * V[age*7+4] * dt
            
            daily_case[age*7+0] = -exposed
            daily_case[age*7+1] = exposed
            daily_case[age*7+2] = as_infection
            daily_case[age*7+3] = m_infection
            daily_case[age*7+4] = mtos_infection
            daily_case[age*7+5] = as_recovery + m_recovery + s_recovery
            daily_case[age*7+6] = m_death + s_death

            Y[age*7+0] = max(V[age*7+0] - exposed, 0)
            Y[age*7+1] = V[age*7+1] + exposed - infection
            Y[age*7+2] = V[age*7+2] + as_infection - as_recovery
            Y[age*7+3] = V[age*7+3] + m_infection - m_recovery - mtos_infection - m_death
            Y[age*7+4] = V[age*7+4] + mtos_infection - s_recovery - s_death
            Y[age*7+5] = V[age*7+5] + as_recovery + m_recovery + s_recovery
            Y[age*7+6] = V[age*7+6] + m_death + s_death
        return Y, daily_case
    
    def split_result(result, idx, u):
        for age in range(n_age):
            result[0][idx, age] = u[age*7+0]
            result[1][idx, age] = u[age*7+1]
            result[2][idx, age] = u[age*7+2]
            result[3][idx, age] = u[age*7+3]
            result[4][idx, age] = u[age*7+4]
            result[5][idx, age] = u[age*7+5]
            result[6][idx, age] = u[age*7+6]
        return result

    u = initial_state
    dt = 1 
    duration = duration - 1
    tl = int(duration*(1/dt))+1
    t = np.linspace(0, duration, tl)
    S = np.zeros((tl, n_age))
    E = np.zeros((tl, n_age))
    Ias = np.zeros((tl, n_age))
    Im = np.zeros((tl, n_age))
    Is = np.zeros((tl, n_age))
    R = np.zeros((tl, n_age))
    D = np.zeros((tl, n_age))
    result = [S, E, Ias, Im, Is, R, D]
    result = split_result(result, 0, u)
         
    d_s = np.zeros((tl, n_age))
    d_e = np.zeros((tl, n_age))
    d_ias = np.zeros((tl, n_age))
    d_im = np.zeros((tl, n_age))
    d_is = np.zeros((tl, n_age))
    d_r = np.zeros((tl, n_age))
    d_d = np.zeros((tl, n_age)) 
    daily_all = [d_s, d_e, d_ias, d_im, d_is, d_r, d_d]
    
    control_time = []

    for j in range(1, tl):
        u, daily = seirmicr(u, dt, t[j], daily_all, control_time)
        result = split_result(result, j, u)
        daily_all = split_result(daily_all, j, daily)

    return result, daily_all, control_time
