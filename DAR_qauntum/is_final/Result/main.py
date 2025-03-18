import scipy.io as scio
import numpy as np
from analyse import calculate_energy
import os

#dataFile1='00306147-VQECImaginaryEvolutionSIQ_ZC-20230106115538_RawData.mat'
#dataFile2='00306149-VQECImaginaryEvolutionSIQ_ZC-20230106125345_RawData.mat'
#dataFile3='00306150-VQECImaginaryEvolutionSIQ_ZC-20230106125611_RawData.mat'
#dataFile4='00306151-VQECImaginaryEvolutionSIQ_ZC-20230106125834_RawData.mat'
#dataFile5='00306152-VQECImaginaryEvolutionSIQ_ZC-20230106130100_RawData.mat'

filename='energy.txt'
with open(filename,'a') as f:
    np.savetxt(f,[['# index    orig_data      orig+rem     orig+rem+cf']],fmt='%s')

#data=scio.loadmat(dataFile1)
#calculate_energy(data,filename)

path=os.getcwd()
f_list=os.listdir(path)

for i in f_list:
    if os.path.splitext(i)[1]=='.mat':
        print('path: ',i)
        data=scio.loadmat(i)
        calculate_energy(data,filename)
