# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 19:49:00 2020
This code is used to implement SIR model with contact matrices.
@author: Jason
"""
import datetime
import numpy as np
import pickle
import os

from function_SEImIcR import model
from class_covid import COVID19
from class_region import Region
from function_postprocessing import sumCaseByAge, paraPrint

"""Fact parameters and model parameters from papers."""
covid = COVID19()
EpidemicParameters = covid.getEpidemicParameters()
paraPrint(EpidemicParameters)

"""Region settings."""
nyPopPath = './data/ny_population.csv' 
nyCMPath = './data/CM/'
nySafeGraphIndexPath = './data/ny_index.csv'
nyTotalInitialInfections = 1000

"""Save path."""
savePath = './simResultsPkl/'
if not os.path.exists(savePath):
    os.mkdir(savePath)

"""Control infomation"""
controlInfoReal = [['2020-01-01', '2020-03-07', {'school': 0.8, 'social distancing': 0.8, 'work': 0.8}],
                ['2020-03-07', '2020-03-17', {'school': 0.8, 'social distancing': 0.5, 'work': 0.5}],
                ['2020-03-17', '2020-03-24', {'school': 0, 'social distancing': 0.4, 'work': 0.5}],
                ['2020-03-24', '2020-06-08', {'school': 0, 'social distancing': 0.3, 'home': 1, 'work': 0.3}],
                ['2020-06-08', '2020-06-22', {'school': 0, 'social distancing': 0.325, 'home': 1, 'work': 0.325}],
                ['2020-06-22', '2020-07-06', {'school': 0, 'social distancing': 0.35, 'home': 1, 'work': 0.35}],
                ['2020-07-06', '2020-07-20', {'school': 0, 'social distancing': 0.375, 'home': 1, 'work': 0.375}],
                ['2020-07-20', '2020-09-01', {'school': 0, 'social distancing': 0.4, 'home': 1, 'work': 0.4}],
                ['2020-09-01', '2020-12-31', {'school': 0.5, 'social distancing': 0.5, 'home': 1, 'work': 0.5}]]
nyReal = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoReal, useSafeGraphIndex=False)
nyRealParameters = nyReal.getRegionParameters()
nyRealCumCase, nyRealDailyCase = model(EpidemicParameters, nyRealParameters)
pickle.dump([nyRealCumCase, nyRealDailyCase], open(savePath+'nyReal_fourPhases.pkl', 'wb'))

controlInfoNothing = [['2020-01-01', '2020-12-31', {'school': 1, 'social distancing': 1}]]
nyNothing = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoNothing, useSafeGraphIndex=False)
nyNothingParameters = nyNothing.getRegionParameters()
nyNothingCumCase, nyNothingDailyCase = model(EpidemicParameters, nyNothingParameters)
pickle.dump([nyNothingCumCase, nyNothingDailyCase], open(savePath+'nyNothing.pkl', 'wb'))

controlInfoSchool = [['2020-01-01', '2020-12-31', {'school': 0, 'social distancing': 1}]]
nySchool = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoSchool, useSafeGraphIndex=False)
nySchoolParameters = nySchool.getRegionParameters()
nySchoolCumCase, nySchoolDailyCase = model(EpidemicParameters, nySchoolParameters)
pickle.dump([nySchoolCumCase, nySchoolDailyCase], open(savePath+'nySchool.pkl', 'wb'))

controlInfoSocialDistance = [['2020-01-01', '2020-12-31', {'school': 0.5, 'social distancing': 0.5, 'work': 0.5}]]
nySocialDistance = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoSocialDistance, useSafeGraphIndex=False)
nySocialDistanceParameters = nySocialDistance.getRegionParameters()
nySocialDistanceCumCase, nySocialDistanceDailyCase = model(EpidemicParameters, nySocialDistanceParameters)
pickle.dump([nySocialDistanceCumCase, nySocialDistanceDailyCase], open(savePath+'nySocialDistance.pkl', 'wb'))

controlInfoElderly = [['2020-01-01', '2020-12-31', {'school': 1, 'social distancing': 1, 'elder': 0.5}]]
nyElderly = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoElderly, useSafeGraphIndex=False)
nyElderlyParameters = nyElderly.getRegionParameters()
nyElderlyCumCase, nyElderlyDailyCase = model(EpidemicParameters, nyElderlyParameters)
pickle.dump([nyElderlyCumCase, nyElderlyDailyCase], open(savePath+'nyElderly.pkl', 'wb'))

controlInfoRelaxJuly = [['2020-01-01', '2020-03-07', {'school': 0.8, 'social distancing': 0.8, 'work': 0.8}],
                ['2020-03-07', '2020-03-17', {'school': 0.8, 'social distancing': 0.5, 'work': 0.5}],
                ['2020-03-17', '2020-03-24', {'school': 0, 'social distancing': 0.4, 'work': 0.5}],
                ['2020-03-24', '2020-06-08', {'school': 0, 'social distancing': 0.3, 'home': 1, 'work': 0.3}],
                ['2020-06-08', '2020-06-22', {'school': 0, 'social distancing': 0.325, 'home': 1, 'work': 0.325}],
                ['2020-06-22', '2020-07-01', {'school': 0, 'social distancing': 0.35, 'home': 1, 'work': 0.35}],
                ['2020-07-01', '2020-09-01', {'school': 0, 'social distancing': 0.5, 'home': 1, 'work': 0.65}],
                ['2020-09-01', '2020-12-31', {'school': 0.5, 'social distancing': 0.5, 'home': 1, 'work': 0.5}]]
nyRelaxJuly = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoRelaxJuly, useSafeGraphIndex=False)
nyRelaxJulyParameters = nyRelaxJuly.getRegionParameters()
nyRelaxJulyCumCase, nyRelaxJulyDailyCase = model(EpidemicParameters, nyRelaxJulyParameters)
pickle.dump([nyRelaxJulyCumCase, nyRelaxJulyDailyCase], open(savePath+'nyRelaxJuly.pkl', 'wb'))

controlInfoRelaxAugust = [['2020-01-01', '2020-03-07', {'school': 0.8, 'social distancing': 0.8, 'work': 0.8}],
                ['2020-03-07', '2020-03-17', {'school': 0.8, 'social distancing': 0.5, 'work': 0.5}],
                ['2020-03-17', '2020-03-24', {'school': 0, 'social distancing': 0.4, 'work': 0.5}],
                ['2020-03-24', '2020-06-08', {'school': 0, 'social distancing': 0.3, 'home': 1, 'work': 0.3}],
                ['2020-06-08', '2020-06-22', {'school': 0, 'social distancing': 0.325, 'home': 1, 'work': 0.325}],
                ['2020-06-22', '2020-07-06', {'school': 0, 'social distancing': 0.35, 'home': 1, 'work': 0.35}],
                ['2020-07-06', '2020-07-20', {'school': 0, 'social distancing': 0.375, 'home': 1, 'work': 0.375}],
                ['2020-07-20', '2020-08-01', {'school': 0, 'social distancing': 0.4, 'home': 1, 'work': 0.4}],
                ['2020-08-01', '2020-12-31', {'school': 0.5, 'social distancing': 0.5, 'home': 1, 'work': 0.5}]]
nyRelaxAugust = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoRelaxAugust, useSafeGraphIndex=False)
nyRelaxAugustParameters = nyRelaxAugust.getRegionParameters()
nyRelaxAugustCumCase, nyRelaxAugustDailyCase = model(EpidemicParameters, nyRelaxAugustParameters)
pickle.dump([nyRelaxAugustCumCase, nyRelaxAugustDailyCase], open(savePath+'nyRelaxAugust.pkl', 'wb'))

controlInfoRelaxOct = [['2020-01-01', '2020-03-07', {'school': 0.8, 'social distancing': 0.8, 'work': 0.8}],
                ['2020-03-07', '2020-03-17', {'school': 0.8, 'social distancing': 0.5, 'work': 0.5}],
                ['2020-03-17', '2020-03-24', {'school': 0, 'social distancing': 0.4, 'work': 0.5}],
                ['2020-03-24', '2020-06-08', {'school': 0, 'social distancing': 0.3, 'home': 1, 'work': 0.3}],
                ['2020-06-08', '2020-06-22', {'school': 0, 'social distancing': 0.325, 'home': 1, 'work': 0.325}],
                ['2020-06-22', '2020-07-06', {'school': 0, 'social distancing': 0.35, 'home': 1, 'work': 0.35}],
                ['2020-07-06', '2020-07-20', {'school': 0, 'social distancing': 0.375, 'home': 1, 'work': 0.375}],
                ['2020-07-20', '2020-09-01', {'school': 0, 'social distancing': 0.4, 'home': 1, 'work': 0.4}],
                ['2020-09-01', '2020-10-01', {'school': 0, 'social distancing': 0.4, 'home': 1, 'work': 0.4}],
                ['2020-10-01', '2020-12-31', {'school': 0.5, 'social distancing': 0.5, 'home': 1, 'work': 0.5}]]
nyRelaxOct = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, nySafeGraphIndexPath, controlInfoRelaxOct, useSafeGraphIndex=False)
nyRelaxOctParameters = nyRelaxOct.getRegionParameters()
nyRelaxOctCumCase, nyRelaxOctDailyCase = model(EpidemicParameters, nyRelaxOctParameters)
pickle.dump([nyRelaxOctCumCase, nyRelaxOctDailyCase], open(savePath+'nyRelaxOct.pkl', 'wb'))



    