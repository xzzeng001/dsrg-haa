import numpy as np
import pyscf
from pyscf import gto,lo,mp
from pyscf import scf,cc
from pyscf import mcscf,mrpt
from pyscf import fci
from pyscf.mcscf import avas
import scipy
import sys
import psutil
import copy

from aug_hessian import orbital_rotation_ah
from gen_ham import generate_ham
from mbe_fci import mbe_fci_corr

#=============FERROCENE====================

mol = gto.Mole()
mol.atom = 'n2.xyz'
mol.basis = 'cc-pvtz'
mol.spin = 0
mol.build()
mol.verbose = 4

mf = scf.RHF(mol)
mf.kernel()
e_hf=mf.e_tot

mo_energy=mf.mo_energy

#mymp = mp.MP2(mf).run()
#mycc=cc.CCSD(mf).run()
#noons, mo_coeff1 = mcscf.addons.make_natural_orbitals(mymp)

#mo=mf.mo_coeff=mo_coeff1

from pyscf import lo
mo_coeff1=mf.mo_coeff
# localization the orbial
new_coeff = lo.Boys(mol, mo_coeff1[:,:int(mol.nelectron/2)]).kernel()
new_coeff2 = lo.Boys(mol, mo_coeff1[:,int(mol.nelectron/2):]).kernel()
mo_coeff2 = np.hstack([new_coeff,new_coeff2])
mf.mo_coeff=mo_coeff2

mo=mf.mo_coeff

from mokit.lib.py2fch_direct import fchk
fchk(mf, 'n2_cmo_boys.fch',density=True)

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
#mycas.max_cycle = 500
mycas.frozen=frozen_lst
orb_indice=[i+1 for i in caslst]
mo1 = mycas.sort_mo(orb_indice)
e_tot, e_cas, fcivec, mo_coeff, mo_energy=mycas.kernel(mo1)
base_e_corr=e_tot-e_hf

mbe_fci_corr(mf,e_tot,base_e_corr,ind_homo,norb,nelec,caslst,vir_index,frozen_lst,mo_coeff,2)

import sys
sys.exit(0)

all_active_lst=[]
ind1,ee,ind1_corr=np.loadtxt('one-body.txt',delimiter=",",unpack=True)
ind2,ind3,ee_2,ind2_corr=np.loadtxt('two-body.txt',delimiter=",",unpack=True)
tol=0.002

all_active_lst=copy.copy(caslst)
for ii,rr in enumerate(ind1_corr):
    if abs(rr) > tol and int(ind1[ii]) not in all_active_lst:
        all_active_lst.append(int(ind1[ii]))

for ii,rr in enumerate(ind2_corr):
    if abs(rr) > tol:
        if int(ind2[ii]) not in all_active_lst:
            all_active_lst.append(int(ind2[ii]))
        if int(ind3[ii]) not in caslst:
            all_active_lst.append(int(ind3[ii]))

n_all_active_orb=len(all_active_lst)

n_all_active_elec=0
for ii in all_active_lst:
    if ii < ind_homo:
        n_all_active_elec += 2

#######################################
# generate hamilonian
#######################################
mycas1 = mcscf.CASCI(mf,n_all_active_orb,n_all_active_elec)
mycas1.max_cycle = 500
mycas1.frozen=frozen_lst

orb_indice=[i+1 for i in all_active_lst]
mo1 = mycas1.sort_mo(orb_indice,mo_coeff=mo_coeff)
e_tot, e_cas, fcivec, mo_coeff, mo_energy=mycas1.kernel(mo1)

generate_ham(mycas1,mo_coeff,mycas1.ncas)

print('energy CASCI: ',e_tot)
base_e_corr=e_tot-e_hf
print('base correlation energy :', base_e_corr)
sys.stdout.flush()


'''
#######################################
# orbital roation
#######################################
ci0=None
eris = mycas.ao2mo(mo_coeff)      
e_tot, e_cas, fcivec = mycas.casci(mo_coeff, ci0, eris)
casdm1, casdm2 = mycas.fcisolver.make_rdm12(fcivec, mycas.ncas, mycas.nelecas)
mo_coeff=orbital_rotation_ah(mycas,mf.mo_coeff,fcivec,eris,casdm1,casdm2,e_cas)
e_tot, e_cas, fcivec=mycas.casci(mo_coeff)

generate_ham(mycas,mo_coeff,mycas.ncas)

print('energy CASCI: ',e_tot)
base_e_corr=e_tot-e_hf

print('base correlation energy :', base_e_corr)
sys.stdout.flush()

########################################
# The NEVPT2 for left orbitals
########################################
ci_nevpt_e1 = mrpt.NEVPT(mycas).kernel()
print('nevpt2: ',ci_nevpt_e1)
sys.stdout.flush()
'''
