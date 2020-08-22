"""
Generate a class represents a specific region (country or city).
@author: Jason 
"""
import pandas as pd
import numpy as np
import datetime

class Region(object):
    def __init__(self, totalInitialInfections, popPath, CMPath, indexPath, controlInfo, useSafeGraphIndex=False):
        self.numAge = 17
        self.popPath = popPath
        self.CMPath = CMPath
        self.indexPath = indexPath
        self.controlInfo = controlInfo
        self.totalInitialInfections = totalInitialInfections
        self.useSafeGraphIndex = useSafeGraphIndex
        
        self.pop = self.getPop()
        self.popRate = self.getPopRate()
        self.originalCM, self.normalizeCM = self.readCM()
        self.initialState = self.getInitialState()
        self.SafeGraphHomeIndex, self.SafeGraphWorkIndex = self.getSafeGraphIndex()
        self.stageList = self.getStageList()
        self.CMList = self.getCMList()
        self.duration = self.getDuration()

    def getRegionParameters(self):
        regionParameters = {'population': self.pop,
                            'populationRate': self.popRate,
                            'contactMatrixList': self.CMList,
                            'duration': self.duration,
                            'initialState': self.initialState,
                            'normalizeContactMatrix': self.normalizeCM}
        return regionParameters
    
    def getPop(self):
        popDF = pd.read_csv(self.popPath)
        p = np.array(popDF['population'])
        return p
    
    def getPopRate(self):
        totalPop = sum(self.pop)
        return [i/totalPop for i in self.pop]
    
    def getInitialState(self):
        initialState = np.zeros((self.numAge*7))
        for age in range(self.numAge):
            initialState[age*7+0] = self.pop[age] - self.popRate[age] * self.totalInitialInfections
            initialState[age*7+2] = 0.4*self.popRate[age] * self.totalInitialInfections
            initialState[age*7+3] = 0.6*self.popRate[age] * self.totalInitialInfections
        return initialState
    
    def getDuration(self):
        start = datetime.datetime.strptime(self.controlInfo[0][0], '%Y-%m-%d')
        end = datetime.datetime.strptime(self.controlInfo[-1][1], '%Y-%m-%d')
        duration = (end - start).days
        return duration
    
    def readCM(self):
        def normalize(cm):
            popRate = np.transpose(np.matrix(self.popRate))
            return (cm + np.multiply(np.transpose(cm), np.dot(popRate, np.transpose(1/popRate))))/2
        def matrixExpand(m):
            m = np.row_stack((m, m[15]))
            m = np.column_stack((m, m[:,15]))
            return m
        home = matrixExpand(np.matrix(pd.read_csv(self.CMPath + 'home.csv')))
        school = matrixExpand(np.matrix(pd.read_csv(self.CMPath + 'school.csv')))
        work = matrixExpand(np.matrix(pd.read_csv(self.CMPath + 'work.csv')))
        other = matrixExpand(np.matrix(pd.read_csv(self.CMPath + 'other.csv')))
        originalCM = [home, school, work, other]      
        home_nor = normalize(home)
        school_nor = normalize(school)
        work_nor = normalize(work)
        other_nor = normalize(other)
        normalizeCM = [home_nor, school_nor, work_nor, other_nor]
        return originalCM, normalizeCM

    def getWeightedCM(self, weightDict):    
        CM = np.dot(weightDict['home'], self.normalizeCM[0]) + np.dot(weightDict['school'], self.normalizeCM[1]) \
            + np.dot(weightDict['work'], self.normalizeCM[2]) + np.dot(weightDict['other'], self.normalizeCM[3])
        return CM 
    
    def getCMList(self):
        CMList = []
        for weightDict in self.stageList:
            CM = self.getWeightedCM(weightDict)
            CMList.append(CM)
        return CMList

    def getSafeGraphIndex(self):
        indexDF = pd.read_csv(self.indexPath)
        SafeGraphHomeIndex = np.array(indexDF['home_change'])
        SafeGraphWorkIndex = np.array(indexDF['work_change'])
        return SafeGraphHomeIndex, SafeGraphWorkIndex               
    
    def getStageList(self):
        """       
        Arguments:
            info {list} -- govenment control information
                info[0]: [start (datetime), end(datetime), control(dict)]
                    control: {'school': 0.5; 'social distancing': 0.6}
        """
        stageList = []
        count = 0
        for idx, stage in enumerate(self.controlInfo):
            weightSchool, weightHome, weightWork, weightOther = 1, 1, 1, 1
            weightElder = 0
            stageStart = datetime.datetime.strptime(stage[0], '%Y-%m-%d')
            stageEnd = datetime.datetime.strptime(stage[1], '%Y-%m-%d')
            SafeGraphIndexStart = datetime.datetime.strptime('2020-02-01', '%Y-%m-%d')
            duration = (stageEnd - stageStart).days
            for key, value in stage[2].items():
                if key == 'school':
                    weightSchool = value
                if key == 'social distancing':
                    weightOther = value
                if key == 'home':
                    weightHome = value
                if key == 'work':
                    weightWork = value
                if key == 'elder':
                    weightElder = value           
            for d in range(duration):
                weightDict = dict()
                if self.useSafeGraphIndex or count > len(self.SafeGraphHomeIndex)-1 or stageEnd <= SafeGraphIndexStart:
                    # use the given weights for home and work
                    weightDict['home'] = np.diag(list(np.repeat(weightHome, self.numAge)))
                    weightDict['work'] = np.diag(list(np.repeat(weightWork, self.numAge)))
                else:
                    # use the index from SafeGraph
                    weightDict['home'] = np.diag(list(np.repeat(self.SafeGraphHomeIndex[count], self.numAge)))
                    weightDict['work'] = np.diag(list(np.repeat(self.SafeGraphWorkIndex[count], self.numAge)))    
                weightDict['school'] = np.diag(list(np.repeat(weightSchool, self.numAge)))
                weightDict['other'] = np.diag(list(np.repeat(weightOther, self.numAge)))
                if weightElder != 0:
                    # the social distancing for the elder is special
                    weightDict['home'] = np.diag(list(np.repeat(1, self.numAge)))
                    weightDict['work'] = np.diag(list(np.repeat(1, 11)) + list(np.repeat(weightElder, 6)))
                    weightDict['school'] = np.diag(list(np.repeat(1, 11)) + list(np.repeat(weightElder, 6)))
                    weightDict['other'] = np.diag(list(np.repeat(1, 11)) + list(np.repeat(weightElder, 6)))
                stageList.append(weightDict)
                count += 1
        return stageList


            
