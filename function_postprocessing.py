# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:54:46 2020

@author: Jason
"""

def sumCaseByAge(cumCase, stateList, ageList, pop_norm=False):
    duration = len(cumCase[stateList[0]])
    sumCumCase = []
    pop_list = [551853, 474248, 482871, 457003, 577721,
                810534, 741552, 629457, 551853,
                551853, 551853, 517362, 491494,
                396644, 293172, 224190, 327662]
    for idx in range(duration):
        sumIdx = 0
        for state in stateList:
            for age in ageList:
                if pop_norm:
                    sumIdx += cumCase[state][idx, age]/pop_list[age]
                else:
                    sumIdx += cumCase[state][idx, age]
        sumCumCase.append(sumIdx)
    return sumCumCase

def paraPrint(epidemicParameters):
    print('COVID Parameters:')
    for key, value in epidemicParameters.items():
        print(key,": ", value)