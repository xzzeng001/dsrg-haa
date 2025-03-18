import scipy.io as scio
import numpy as np
from analyse import calculate_energy
import os

filename='energy.txt'
with open(filename,'a') as f:
    np.savetxt(f,[['# index    orig_data      orig+rem     orig+rem+cf     VQE_no_noise']],fmt='%s')

path=os.getcwd()
f_list=os.listdir(path)

for i in f_list:
    if os.path.splitext(i)[1]=='.mat':
        print('path: ',i)
        data=scio.loadmat(i)
        calculate_energy(data,filename,1)

#datafile='00306076-VQECImaginaryEvolutionSIQ_ZC-20230105201313_RawData.mat'
#data=scio.loadmat(datafile)
#params=data['paramsList']
##print('params-1:',params[0])
#
#energy=main_vqe(1,params[0])
#
##print('energy: ',energy)
#
