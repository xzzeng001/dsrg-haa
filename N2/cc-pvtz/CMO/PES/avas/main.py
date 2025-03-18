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

file_name='energy_pes.txt'
np.savetxt(file_name,['# bond-length  CASCI  CASSCF  CASSCF+NEVPT2'],fmt='%s')
current_path = os.getcwd()

for rr in np.array([0.6,0.8,1.0,1.09,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0]):

    # mkdir the work directory
    file_path0 = current_path+'/'+str(rr.round(2))
    file_mo_energy=file_path0+'/mo_energy.txt'
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

#    # for mp2 energy
#    mymp = mp.MP2(mf).run()
#    noons, mo_coeff1 = mcscf.addons.make_natural_orbitals(mymp)
#    mf.mo_coeff=mo_coeff1

    #from pyscf.mcscf import avas
    ao_labels = ['N 2p']
    avas_obj = avas.AVAS(mf, ao_labels)
    avas_obj.kernel()
        
    nocc=len(avas_obj.occ_weights)
    caslst=[]
    for ii,rtmp in enumerate(avas_obj.occ_weights):
        if rtmp > avas_obj.threshold:
            caslst.append(ii)
        
    for ii,rtmp in enumerate(avas_obj.vir_weights):
        if rtmp > avas_obj.threshold:
            caslst.append(ii+nocc)
    
    n_orb=len(caslst)
        
    with open(file_orb_list,'a') as f:
        np.savetxt(f,[caslst],fmt='%d')
    
    mycas1 = mcscf.CASCI(mf,avas_obj.ncas,avas_obj.nelecas)
    ## init_cas add 1 for indice
    orb_indice=[i+1 for i in caslst]

    with open(file_orb_list,'a') as f:
        np.savetxt(f,[orb_indice],fmt='%d')

    mo1 = mycas1.sort_mo(orb_indice)
    #
    e_tot_1, e_cas, fcivec, mo, mo_energy=mycas1.kernel(mo1)

    mycas2 = mcscf.CASSCF(mf,avas_obj.ncas,avas_obj.nelecas)
    mo2 = mycas2.sort_mo(orb_indice)
    e_tot_2, e_cas, fcivec, mo, mo_energy=mycas2.kernel(mo2)
    #       
        
    ci_nevpt_e1 = mrpt.NEVPT(mycas2).kernel()
 
    with open(file_name,'a') as f:
        np.savetxt(f,[[rr, e_tot_1, e_tot_2, e_tot_2+ci_nevpt_e1]]) 
