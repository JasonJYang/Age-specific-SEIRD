import pickle
import numpy as np
from function_postprocessing import sumCaseByAge

def checkEffect(real, simulate, stateList):
    real = sumCaseByAge(real, stateList, list(range(17)))
    simulate = sumCaseByAge(simulate, stateList, list(range(17)))
    real_sum = sum(real)
    simulate_sum = sum(simulate)
    return simulate_sum <= real_sum
    
def economicEffect(weight, policy):
    return np.dot(np.array(weight).transpose(), np.array(policy))
    
nyRealCumCase, nyRealDailyCase = pickle.load(open('./simResultsPkl/nyReal_fourPhases.pkl', 'rb'))

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

def best_search(weight, policy_list):                
    maxScore = 0
    bestPolicy = []
    for policy in policy_list:
        score = economicEffect(weight, policy)
        if score > maxScore:
            maxScore = score
            bestPolicy = policy[:]
    return bestPolicy
  
weight_identity = [1]*4
weight_heterogeneity = [1, 0, 5, 2.5] 
print('Identity weight')            
print('Best control for infections', best_search(weight_identity, effectTotalList))
print('Best control for dealh', best_search(weight_identity, effectDeathList))
print('Best control for severe infections', best_search(weight_identity, effectSevereList))
print('Heterogeneity weights')
print('Best control for infections', best_search(weight_heterogeneity, effectTotalList))
print('Best control for dealh', best_search(weight_heterogeneity, effectDeathList))
print('Best control for severe infections', best_search(weight_heterogeneity, effectSevereList))