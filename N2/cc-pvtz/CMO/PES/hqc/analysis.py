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

file_name='energy_pes.txt'
np.savetxt(file_name,['# bond-length  CASCI  CASSCF  CASSCF+NEVPT2'],fmt='%s')
current_path = os.getcwd()

for rr in np.array([0.6,0.8,1.0,1.09,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0]):

    # mkdir the work directory
    file_path0 = current_path+'/'+str(rr.round(2))
    file_mo_energy=file_path0+'/mo_energy.txt'
    file_one_body=file_path0+'/one_body.txt'
    file_two_body=file_path0+'/two_body.txt'
    file_orb_list=file_path0+'/orb_list.txt'
#    os.mkdir(file_path0) 

    mol = gto.Mole()
    mol.atom = [["N", [0., 0., 0.0]],["N", [0., 0., rr]]]
    mol.basis = 'cc-pvtz'
    mol.spin = 0
    mol.build()
    mol.verbose = 4

    mf = scf.RHF(mol)
    mf.kernel()
    e_hf=mf.e_tot

    caslst=np.loadtxt(file_orb_list)

    ind_homo=sum(mol.nelec) // 2
    n_orb=mol.nao_nr()

    all_active_lst=[int(i) for i in caslst]
    
    n_orb=len(all_active_lst)
    
    n_elec=0
    for ii in all_active_lst:
        if ii < ind_homo:
            n_elec += 2
    
    mycas1 = mcscf.CASCI(mf,n_orb,n_elec)
    
    orb_indice=[i+1 for i in all_active_lst]

    mo1 = mycas1.sort_mo(orb_indice)
    e_casci, e_cas, fcivec, mo_coeff, mo_energy=mycas1.kernel(mo1)

    mycas2 = mcscf.CASSCF(mf,n_orb,n_elec)
    mo2 = mycas2.sort_mo(orb_indice)
    e_casscf, e_cas, fcivec, mo_coeff, mo_energy=mycas2.kernel(mo2)
    
    ci_nevpt_e1 = mrpt.NEVPT(mycas2).kernel()    
  
    with open(file_name,'a') as f:
        np.savetxt(f,[[rr, e_casci,e_casscf,ci_nevpt_e1+e_casscf]]) 
