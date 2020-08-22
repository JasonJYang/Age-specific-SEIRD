# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:20:05 2020
Possible solution save.

@author: Jason
"""
import numpy as np
import os
import pickle
from function_SEImIcR import model
from class_covid import COVID19
from class_region import Region
from function_postprocessing import paraPrint

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
savePath = './OptResultsPkl/'
if not os.path.exists(savePath):
    os.mkdir(savePath)

weightWorkList = weightSchoolList = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
weightOtherList = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
for weightSchool in weightSchoolList:
        for weightWork in weightWorkList:
            for weightOther in weightOtherList:
                controlInfoSimulate = [['2020-01-01', '2020-12-31', {'school': weightSchool,
                                                                       'home': 1,
                                                                       'work': weightWork,
                                                                       'social distancing': weightOther}]]
                nySimulate = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, 
                                    nySafeGraphIndexPath, controlInfoSimulate, useSafeGraphIndex=False)
                nySimulateParameters = nySimulate.getRegionParameters()
                nySimulateCumCase, nySimulateDailyCase = model(EpidemicParameters, nySimulateParameters)
                record = 'school_'+str(weightSchool)+'_'+'work_'+str(weightWork)+'_'+'other_'+str(weightOther)
                savePathCur = savePath + record
                output = open(savePathCur+'.pkl', 'wb')
                pickle.dump([nySimulateCumCase, nySimulateDailyCase], output)
                output.close()
