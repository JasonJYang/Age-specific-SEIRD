"""
COVID information
"""
import math
import numpy as np

class COVID19(object):
    def __init__(self):
        self.R0 = 2.0
        self.numAge = 17
        self.latentPeriod = 6
        self.infectiousPeriod = 4.6
        self.mildPeriod = 10
        self.onsetToSeverePeriod = 3
        self.onsetToDeathPeriod = 10
        self.onsetToDischargePeriod = 16.25       
        self.confirmRateList = np.repeat(0.4, self.numAge)
        self.severeRateList = [0.025, 0.025, 0.025, 0.025,
                               0.32, 0.32, 0.32, 0.32,
                               0.32, 0.32, 0.32, 0.32,
                               0.64, 0.64, 0.64, 0.64,
                               0.64]
        self.crudeMortalityRateList = [0.0215, 0.0215, 0.0215, 0.0782,
                                       0.0782, 0.0782, 0.0782, 0.0782,
                                       0.0782, 0.2257, 0.2257, 0.2257,
                                       0.2257, 0.46, 0.46, 0.6085,
                                       0.6085]

    def getEpidemicParameters(self):        
        self.transmissionRateList = np.repeat(self.R0*(1/self.infectiousPeriod), self.numAge)        
        self.infectiousRateList = np.repeat(1/self.latentPeriod, self.numAge) # alpha     
        self.mildToSevereRateList = [i*(1/self.onsetToSeverePeriod) for i in self.severeRateList]        
        self.asymRecoverRateList = np.repeat(1/self.infectiousPeriod, self.numAge)              
        self.mildRecoverRateList = np.repeat(1/self.mildPeriod, self.numAge) # gamma       
        self.severeRecoveRateList = np.repeat(1/(self.onsetToDischargePeriod-self.onsetToSeverePeriod), self.numAge) # gamma      
        self.mortalityRateList = [i*(1/(self.onsetToDeathPeriod-self.onsetToSeverePeriod)) for i in self.crudeMortalityRateList]
        EpidemicParameters = {'R0': self.R0,
                          'transmissionRateList': self.transmissionRateList,
                          'infectiousRateList': self.infectiousRateList,
                          'mildToSevereRateList': self.mildToSevereRateList,
                          'asymRecoverRateList': self.asymRecoverRateList,
                          'mildRecoverRateList': self.mildRecoverRateList,
                          'severeRecoveRateList': self.severeRecoveRateList,
                          'mortalityRateList': self.mortalityRateList,
                          'confirmRateList': self.confirmRateList}
        return EpidemicParameters

    def setR0(self, R0):
        self.R0 = R0

    def setLatentPeriod(self, latentPeriod):
        self.latentPeriod = latentPeriod

    def setInfectiousPeriod(self, infectiousPeriod):
        self.infectiousPeriod = infectiousPeriod
    
    def setMildPeriod(self, mildPeriod):
        self.mildPeriod = mildPeriod
    
    def setOnsetToSeverePeriod(self, onsetToSeverePeriod):
        self.onsetToSeverePeriod = onsetToSeverePeriod
    
    def setOnsetToDeathPeriod(self, onsetToDeathPeriod):
        self.onsetToDeathPeriod = onsetToDeathPeriod
    
    def setOnsetToDischargePeriod(self, onsetToDischargePeriod):
        self.onsetToDischargePeriod = onsetToDischargePeriod
    
    def setConfirmRateList(self, confirmRateList):
        self.confirmRateList = confirmRateList
    
    def setSevereRateList(self, severeRateList):
        self.severeRateList = severeRateList
    
    def setCrudeMortalityRateList(self, crudeMortalityRateList):
        self.crudeMortalityRateList = crudeMortalityRateList
    

        
        