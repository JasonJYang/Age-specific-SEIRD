import pickle
import numpy as np
from function_postprocessing import sumCaseByAge
from class_region import Region

def checkEffect(real, simulate, stateList):
    real = sumCaseByAge(real, stateList, list(range(17)))
    simulate = sumCaseByAge(simulate, stateList, list(range(17)))
    real_sum = sum(real)
    simulate_sum = sum(simulate)
    return simulate_sum <= real_sum
    
def economicEffect(weight, policy, CM):
    contact = np.array(policy) * CM
    return np.dot(np.array(weight).transpose(), contact)

def best_search(weight, policy_list, CM):                
    maxScore = 0
    bestPolicy = []
    for policy in policy_list:
        score = economicEffect(weight, policy, CM)
        if score > maxScore:
            maxScore = score
            bestPolicy = policy[:]
    return bestPolicy
    
nyRealCumCase, nyRealDailyCase = pickle.load(open('./simResultsPkl/nyReal_fourPhases.pkl', 'rb'))

"""Region settings."""
nyPopPath = './data/ny_population.csv' 
nyCMPath = './data/CM/'
nySafeGraphIndexPath = './data/ny_index.csv'
nyTotalInitialInfections = 1000
controlInfoSimulate = [['2020-01-01', '2020-12-31', {'school': 1,
                                                     'home': 1,
                                                     'work': 1,
                                                     'social distancing': 1}]]
nySimulate = Region(nyTotalInitialInfections, nyPopPath, nyCMPath, 
                    nySafeGraphIndexPath, controlInfoSimulate, useSafeGraphIndex=False)
CM = nySimulate.normalizeCM
CM_sum = [np.sum(val) for val in CM]

path0 = './OptResultsPkl/'
weightWorkList = weightSchoolList = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
weightOtherList = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

effectTotalList = []
effectDeathList = []
effectSevereList = []
for weightSchool in weightSchoolList:
    for weightWork in weightWorkList:
        for weightOther in weightOtherList:
            record = 'school_'+str(weightSchool)+'_'+'work_'+str(weightWork)+'_'+'other_'+str(weightOther)
            save_path = path0 + record
            fr = open(save_path+'.pkl','rb')
            nySimulateCumCase, nySimulateDailyCase = pickle.load(fr)
            fr.close()
            
            if checkEffect(nyRealDailyCase, nySimulateDailyCase, [2,3,4]):
                effectTotalList.append([weightSchool, 1, weightWork, weightOther])
                
            if checkEffect(nyRealDailyCase, nySimulateDailyCase, [4]):
                effectSevereList.append([weightSchool, 1, weightWork, weightOther])
               
            if checkEffect(nyRealDailyCase, nySimulateDailyCase, [6]):
                effectDeathList.append([weightSchool, 1, weightWork, weightOther])
  
weight_identity = [1]*4
weight_heterogeneity = [0.3, 0.3, 1, 0.5] 
print('Identity weight')            
print('Best control for infections', best_search(weight_identity, effectTotalList, CM_sum))
print('Best control for death', best_search(weight_identity, effectDeathList, CM_sum))
print('Best control for severe infections', best_search(weight_identity, effectSevereList, CM_sum))
print('Intuitive assumptions')
print('Best control for infections', best_search(weight_heterogeneity, effectTotalList, CM_sum))
print('Best control for death', best_search(weight_heterogeneity, effectDeathList, CM_sum))
print('Best control for severe infections', best_search(weight_heterogeneity, effectSevereList, CM_sum))