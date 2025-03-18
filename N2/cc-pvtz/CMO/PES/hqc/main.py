import numpy as np
import pyscf
from pyscf import gto,lo,mp
from pyscf import scf,cc,ci
from pyscf import mcscf,mrpt
from pyscf import fci
from pyscf.mcscf import avas
import scipy
import sys
import psutil
import copy
import os 

from aug_hessian import orbital_rotation_ah
from gen_ham import generate_ham
from mbe_fci import mbe_fci_corr

#file_name='energy_pes.txt'
#np.savetxt(file_name,['# bond-length   RHF   MP2   CCSD   CISD   CASCI CASSCF  CASSCF+NEVPT2'],fmt='%s')
current_path = os.getcwd()

for rr in np.array([0.6,0.8,1.0,1.09,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0]):

    # mkdir the work directory
    file_path0 = current_path+'/'+str(rr.round(2))
    file_mo_energy=file_path0+'/mo_energy.txt'
    file_one_body=file_path0+'/one_body.txt'
    file_two_body=file_path0+'/two_body.txt'
    file_orb_list=file_path0+'/orb_list.txt'
    os.mkdir(file_path0) 

    mol = gto.Mole()
    mol.atom = [["N", [0., 0., 0.0]],["N", [0., 0., rr]]]
    mol.basis = 'cc-pvtz'
    mol.spin = 0
    mol.build()
    mol.verbose = 4

    mf = scf.RHF(mol)
    mf.kernel()
    e_hf=mf.e_tot

    np.savetxt(file_mo_energy,mf.mo_energy)

    # for mp2 energy
#    mymp = mp.MP2(mf).run()

    # calculate the no orbitals
#    noons, mo_coeff1 = mcscf.addons.make_natural_orbitals(mymp)
#    mf.mo_coeff=mo_coeff1

    #######################################
    # orbital selection based on the MBE
    #######################################
    ind_homo=sum(mol.nelec) // 2
    n_orb=mol.nao_nr()
    all_orb_list=[i for i in range(n_orb)]
    norb=2
    nelec=2
    caslst=[ind_homo-1,ind_homo]
    frozen_lst=[0,1,50,51,52,53,54,55,56,57,58,59]
    vir_index=[]
    for ii in all_orb_list:
        if ii not in caslst and ii not in frozen_lst:
            vir_index.append(ii)
    
    mycas = mcscf.CASCI(mf,norb,nelec)
    mycas.frozen=frozen_lst
    orb_indice=[i+1 for i in caslst]
    mo1 = mycas.sort_mo(orb_indice)
    e_tot, e_cas, fcivec, mo_coeff, mo_energy=mycas.kernel(mo1)
    base_e_corr=e_tot-e_hf
    
    mbe_fci_corr(mf,e_tot,base_e_corr,ind_homo,norb,nelec,caslst,vir_index,frozen_lst,mo_coeff,2,file_one_body,file_two_body)
    
#    all_active_lst=[]
#    ind1,ee,ind1_corr=np.loadtxt(file_one_body,delimiter=",",unpack=True)
#    ind2,ind3,ee_2,ind2_corr=np.loadtxt(file_two_body,delimiter=",",unpack=True)
#    tol=0.005
#    
#    all_active_lst=copy.copy(caslst)
#    for ii,rtmp in enumerate(ind1_corr):
#        if abs(rtmp) > tol and int(ind1[ii]) not in all_active_lst:
#            all_active_lst.append(int(ind1[ii]))
#    
#    for ii,rtmp in enumerate(ind2_corr):
#        if abs(rtmp) > tol:
#            if int(ind2[ii]) not in all_active_lst:
#                all_active_lst.append(int(ind2[ii]))
#            if int(ind3[ii]) not in caslst:
#                all_active_lst.append(int(ind3[ii]))
#    
#    n_all_active_orb=len(all_active_lst)
#    
#    n_all_active_elec=0
#    for ii in all_active_lst:
#        if ii < ind_homo:
#            n_all_active_elec += 2
#    
#    #######################################
#    # generate hamilonian
#    #######################################
#    mycas1 = mcscf.CASCI(mf,n_all_active_orb,n_all_active_elec)
#    mycas1.frozen=frozen_lst
#    
#    orb_indice=[i+1 for i in all_active_lst]
#
#    with open(file_orb_list,'a') as f:
#        np.savetxt(f,[[n_all_active_orb,n_all_active_elec]],fmt='%d,%d')
#
#    with open(file_orb_list,'a') as f:
#        np.savetxt(f,[orb_indice],fmt='%d')
#
#    mo1 = mycas1.sort_mo(orb_indice)
#    e_tot_1, e_cas, fcivec, mo_coeff, mo_energy=mycas1.kernel(mo1)
#
#    mycas2 = mcscf.CASSCF(mf,n_all_active_orb,n_all_active_elec)
#    mo2 = mycas2.sort_mo(orb_indice)
#    e_tot_2, e_cas, fcivec, mo_coeff, mo_energy=mycas2.kernel(mo2)
#    
#    ci_nevpt_e1 = mrpt.NEVPT(mycas2).kernel()    
##    generate_ham(mycas1,mo_coeff,mycas1.ncas)
   
##    with open(file_name,'a') as f:
#        np.savetxt(f,[[rr, e_hf, mp_energy, cc_energy, ci_energy]]) 
