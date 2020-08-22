# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 19:49:00 2020
This code is used to implement SIR model with contact matrices.
@author: Jason
"""
import datetime
import pickle
from collections import defaultdict
import numpy as np
import os
from function_SEImIcR_adaptive import model
from class_covid import COVID19
from class_region import Region
from function_postprocessing import sumCaseByAge, paraPrint

def control_list(control_time):
    result_date = []
    start = current = 0
    while current != len(control_time):
        if current == len(control_time)-1:
            date_tuple = (control_time[start], control_time[current])
            result_date.append(date_tuple) 
        if control_time[current] - control_time[start] != current - start:
            date_tuple = (control_time[start], control_time[current-1])
            result_date.append(date_tuple)
            start = current
        else:
            current += 1
    return result_date

"""Fact parameters and model parameters from papers."""
covid = COVID19()
EpidemicParameters = covid.getEpidemicParameters()
paraPrint(EpidemicParameters)

"""Save path."""
savePath = './simResultsPkl/'
if not os.path.exists(savePath):
    os.mkdir(savePath)

"""Region settings."""
nyPopPath = './data/ny_population.csv' 
nyCMPath = './data/CM/'
nySafeGraphIndexPath = './data/ny_index.csv'
nyTotalInitialInfections = 1000

"""Model"""
controlInfoAdaptive = [['2020-01-01', '2020-12-03', {'school': 1}]]
nyAdaptive = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoAdaptive, useSafeGraphIndex=False)
nyAdaptiveParameters = nyAdaptive.getRegionParameters()

controlSevere100 = {'severe': 100}
nyAdaptiveSevere100CumCase, nyAdaptiveSevere100DailyCase, nyAdaptiveSevere100ControlTime = model(EpidemicParameters,
                                                                                               nyAdaptiveParameters,
                                                                                               controlSevere100,
                                                                                               nyAdaptive)
severe100ControlDate = control_list(nyAdaptiveSevere100ControlTime)
pickle.dump([nyAdaptiveSevere100CumCase, 
             nyAdaptiveSevere100DailyCase, 
             severe100ControlDate], 
             open(savePath+'nyAdaptiveSevere100.pkl', 'wb'))

controlSevere150 = {'severe': 150}
nyAdaptiveSevere150CumCase, nyAdaptiveSevere150DailyCase, nyAdaptiveSevere150ControlTime = model(EpidemicParameters,
                                                                                               nyAdaptiveParameters,
                                                                                               controlSevere150,
                                                                                               nyAdaptive)
severe150ControlDate = control_list(nyAdaptiveSevere150ControlTime)
pickle.dump([nyAdaptiveSevere150CumCase, 
             nyAdaptiveSevere150DailyCase, 
             severe150ControlDate], 
             open(savePath+'nyAdaptiveSevere150.pkl', 'wb'))







    
    
    
    