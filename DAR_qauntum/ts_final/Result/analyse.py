import numpy as np
from qc_fun import _readout_error_mitigation,get_energy_category,get_expectation_value_single,bin2numlist  

def calculate_energy(data,filename):
    ## calculate the energy from original data

    observables_tmp=data['observableH']
    observables=[]
    for rtmp in observables_tmp[0]:
        rr=str(rtmp).split('\'')[1]
        observables.append(rr)
    
    #print('observables: ',observables)
    observableCommute=data['observableHCommuteBasis']
    measurementBasis=data['measurementBasisH']
    bitCount=data['bitCountList']
    REMMatrix=data['REMMatrix']
    const=data['observableProbabilityH']
    cliffordFitParameters=data['cliffordFitParametersH']
    
    ## for energy with REM and CF
    nn=len(bitCount)
    for ii in range(nn):
        res_noise=bitCount[ii]
        Probability=get_expectation_value_single(observables, observableCommute, measurementBasis, res_noise)
        energy_without_rem=np.dot(Probability,const.transpose())
    
        res_noise_with_rem,bitStringCountTmp=_readout_error_mitigation(REMMatrix, res_noise)
        Probability_with_rem=get_expectation_value_single(observables, observableCommute, measurementBasis, res_noise_with_rem)
        energy_with_rem=np.dot(Probability_with_rem,const.transpose())
    
        keyValue = np.array(Probability_with_rem) * np.array(cliffordFitParameters[0]) + np.array(cliffordFitParameters[1])
        keyValue = np.clip(keyValue,-1,1)
        ### 计算能量值
        res_H = np.dot(keyValue,const.transpose())
    
        with open(filename,'a') as f:
            np.savetxt(f,[[ii+1,energy_without_rem,energy_with_rem,res_H]],fmt='%d     %15.8f     %15.8f     %15.8f')
    
