import numpy as np
import pyscf
from pyscf import gto,lo,mp
from pyscf import scf,cc
from pyscf import mcscf,mrpt
from pyscf import fci
import scipy
import sys
import psutil
import copy

from mbe_fci import mbe_fci_corr

#=============FERROCENE====================

mol = gto.Mole()
mol.atom = 'feco.xyz'
mol.basis = {'Fe':'lanl2dz','C':'def2-svp','O':'def2-svp'}
mol.ecp = {'Fe':'lanl2dz'}
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

mo=mf.mo_coeff

#from mokit.lib.py2fch_direct import fchk
#fchk(mf, 'n2_cmo.fch',density=True)

#######################################
# orbital selection based on the MBE
#######################################
ind_homo=sum(mol.nelec) // 2
n_orb=mol.nao_nr()
all_orb_list=[i for i in range(n_orb)]
norb=2
nelec=2
caslst=[ind_homo-1,ind_homo]
frozen_lst=None
vir_index=[]
for ii in all_orb_list:
    if ii not in caslst: #and ii not in frozen_lst:
        vir_index.append(ii)

mycas = mcscf.CASCI(mf,norb,nelec)
#mycas.max_cycle = 500
mycas.frozen=frozen_lst
orb_indice=[i+1 for i in caslst]
mo1 = mycas.sort_mo(orb_indice)
e_tot, e_cas, fcivec, mo_coeff, mo_energy=mycas.kernel(mo1)
base_e_corr=e_tot-e_hf

import time
start_time=time.time()
mbe_fci_corr(mf,e_tot,base_e_corr,ind_homo,norb,nelec,caslst,vir_index,frozen_lst,mo_coeff,2)
end_time=time.time()

print('Total time for MBECAS: ',end_time-start_time)

