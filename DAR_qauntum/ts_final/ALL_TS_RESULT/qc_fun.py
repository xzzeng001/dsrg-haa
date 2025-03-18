import numpy as np
import math

def _readout_error_mitigation(REMMatrix, bitString):#,num_sots):
    '''
    从原始bitstring期望值得到经过读取校准后的bitstring期望值
    REMMatrix：读取矩阵，数据文件中的REMMatrix
    bitString：原始bitstring期望值，数据文件中的bitCountList
    '''
#    qubitsNumber = len(bitString)
    qubitsNumber = int(math.log2(bitString.shape[1]))
#    print('bitString.shape: ',bitString.shape)
#    print('shape[2]: ',bitString.shape[1])

    REMMatrix_kron = [[REMMatrix[0][0],REMMatrix[0][1]],[REMMatrix[1][0],REMMatrix[1][1]]]
    for ii in range(1,qubitsNumber):
        RMatrix_Tmp = [[REMMatrix[ii*2][0],REMMatrix[ii*2][1]],[REMMatrix[ii*2+1][0],REMMatrix[ii*2+1][1]]]
        REMMatrix_kron = np.kron(REMMatrix_kron,RMatrix_Tmp)
    
#    bin2decList = []
#    for ii in range(qubitsNumber):
#        bin2decList.append(2**(qubitsNumber-ii-1))
#        
#    bitStringCountTmp = np.zeros(2**qubitsNumber)
##    num_shots = len(bitString[0])
#    bitStringVector = np.dot(bin2decList,np.array(bitString))
#    for bitStringDec in bitStringVector:
#        bitStringCountTmp[bitStringDec] +=1
#    bitStringCountTmp = np.array(bitStringCountTmp)/num_shots
    bitStringCountTmp=bitString.transpose()
    bitProbability = np.dot(REMMatrix_kron,bitStringCountTmp)
        
    return bitProbability.transpose(),bitStringCountTmp


def get_energy_category(observables, observableCommuteBasisMatrix, measurementBasis, cliffordFitParameters, obserProbability, res_noise):
    '''
    计算经过clifford fit的能量值
    '''
    expValueH = get_expectation_value_single(observables, observableCommuteBasisMatrix, measurementBasis, res_noise)##计算各个observable的期望值
    ### Clifford Fit
    keyValue = np.array(expValueH) * np.array(cliffordFitParameters[0]) + np.array(cliffordFitParameters[1])
    keyValue = np.clip(keyValue,-1,1)
    ### 计算能量值
    res_H = np.dot(keyValue,obserProbability)
    del expValueH
    del keyValue
    gc.collect()
    return res_H

def get_expectation_value_single(observables, observableCommute, measurementBasis, res_noise):
    '''
    从经过读取矫正的bitstring期望值结果，得到observable期望值
    observables：observable的list，数据文件中的observableH
    observableCommute：测量basis和observable之间的对易关系，数据文件中的observableHCommuteBasis
    measurementBasis：测量的basis，数据文件中的measurementBasisH
    res_noise：经过读取矫正的bitstring期望值，由bitCountList和REMMatrix计算获得
    '''
    qubitsNumber = len(observables[0])
    binaryLen = '0' + str(qubitsNumber) + 'b'

    observNumList = []
    commuteNormal = []
    observableCommute = np.array(observableCommute)
    for ii,ob in enumerate(observables):
        observNumList.append(obserc2numlist(ob))
        commuteNormal.append(max([1,sum(observableCommute[ii])]))
        
    binNumList = []
    for decValue in range(2**qubitsNumber):
        binaryBit = format(decValue,binaryLen)
#        print('binaryBit: ',binaryBit)
        binNumList.append(bin2numlist(binaryBit))
#        print('bin2numlist(binaryBit): ',bin2numlist(binaryBit))

#    print('binNumList: ',np.array(binNumList).shape)
#    print('observNumList: ',np.array(observNumList).shape)
    Martrix = 1 - 2*np.mod(np.dot(np.array(binNumList), np.array(observNumList).transpose()),2)
#    print('Martrix: ',Martrix.shape)
#    print('res_noise: ',res_noise.shape)
    expTmp = np.dot(res_noise, Martrix)
    # expValueTmp = np.dot(observableCommute,expTmp)
    # exp_value = np.diagonal(expValueTmp)/commuteNormal
    expValueTmp = []
    for ii in range(len(expTmp[0])):
        resTmp = np.dot(observableCommute[ii,:], expTmp[:,ii])
        expValueTmp.append(resTmp)
    exp_value = np.array(expValueTmp)/commuteNormal
    del expValueTmp
    del expTmp
#    gc.collect()
    return exp_value

def obserc2numlist(observ):
    obnumlist = []
    for idxStr,pauliStr in enumerate(observ):
        if pauliStr == 'I':
            obnumlist.append(0)
        else:
            obnumlist.append(1)
    return obnumlist
            
def bin2numlist(binValue):
    binumlist = []
    for idxStr,binStr in enumerate(binValue):
        if binStr == '0':
            binumlist.append(0)
        else:
            binumlist.append(1)
    return binumlist
