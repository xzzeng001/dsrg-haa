import numpy as np
import os

# do some preparation
current_path = os.getcwd()
file_one_body=current_path+'/one_body.txt'
file_two_body=current_path+'/two_body.txt'
file_orb_list=current_path+'/orb_list.txt'
caslst=[6,7]  # index for HOMO,LUMO

# read the content of orbtial correlation 
ind1,ee,ind1_corr=np.loadtxt(file_one_body,delimiter=",",unpack=True)
ind2,ind3,ee_2,ind2_corr=np.loadtxt(file_two_body,delimiter=",",unpack=True)

ind1_new=[int(ii) for ii in ind1]
ind2_new=[int(ii) for ii in ind2]
ind3_new=[int(ii) for ii in ind3]

ind1_corr_abs=abs(ind1_corr)
ind2_corr_abs=abs(ind2_corr)
nn1=ind1_corr_abs.shape

# sum the total abs energy
e_corr_abs=np.append(ind1_corr_abs,ind2_corr_abs)
e_corr_ind=np.argsort(e_corr_abs)[::-1]
e_tot=sum(e_corr_abs)

tol=0.4
e_new=0

for ind_ in e_corr_ind:
    e_new=e_new+e_corr_abs[ind_]
    if e_new/e_tot < tol:
        if ind_ < nn1:
            new_cas=ind1_new[ind_]
            if new_cas not in caslst:
                caslst.append(new_cas)
        else:
            new_cas=ind2_new[int(ind_-nn1)]
            if new_cas not in caslst:
                caslst.append(new_cas)

            new_cas=ind3_new[int(ind_-nn1)]
            if new_cas not in caslst:
                caslst.append(new_cas)
    else:
        break


with open(file_orb_list,'a') as f:
    np.savetxt(f,[np.sort(caslst)],fmt='%d')

#print('caslst: ', caslst)

