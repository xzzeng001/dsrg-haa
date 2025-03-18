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

mymp = mp.MP2(mf).run()
#mycc=cc.CCSD(mf).run()
noons, mo_coeff1 = mcscf.addons.make_natural_orbitals(mymp)

mo=mf.mo_coeff=mo_coeff1

norb=5
nelec=6

#frozen_lst=[0,1,50,51,52,53,54,55,56,57,58,59]
mycas = mcscf.CASSCF(mf,norb,nelec)
#mycas.frozen=frozen_lst
orb_indice=[4,6,7,8,9]
mo1 = mycas.sort_mo(orb_indice)
e_tot, e_cas, fcivec, mo_coeff, mo_energy=mycas.kernel(mo1)

########################################
# The NEVPT2 for left orbitals
########################################
from pyscf import mrpt
e_corr = mrpt.NEVPT(mycas).kernel()
e_tot2 = e_tot + e_corr
print('nevpt2 and energy: ',e_corr,e_tot2)
sys.stdout.flush()
