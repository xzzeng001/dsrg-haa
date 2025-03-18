import functools
from functools import reduce
import numpy as np
import pyscf
from pyscf import gto,lo,mp
from pyscf import scf,cc
from pyscf.cc import ccsd
from pyscf import mcscf,mrpt
from pyscf import fci,ci
from pyscf.mcscf import avas
import scipy
import sys
import psutil
import copy

import openfermion
from openfermion import MolecularData
import os
import cmath

#=============H10_chain====================

file_name='energy_pes.txt'
np.savetxt(file_name,['# bond-length   RHF   MP2   CCSD   CISD   CASCI   CASSCF   CASSCF+NEVPT2'],fmt='%s')

current_path = os.getcwd()

mol = gto.M()
mol.atom = "ts.xyz"
mol.basis ='6-31G(D)'
mol.charge = 0
mol.spin = 0
mol.verbose = 4
mol.build(parse_arg=False)

# for the RHF energy
mf = scf.RHF(mol)
rhf_energy=mf.kernel()

## for the MP2  energy
#mymp = mp.MP2(mf).run()
#mp_energy=mymp.e_tot
#
## for the CCSD energy
#mycc=cc.CCSD(mf).run()
#cc_energy=mycc.e_tot
#
## for the CISD energy
#myci = ci.CISD(mf).run()
#ci_energy=myci.e_tot

# for CASCI energy
norb=4
nelec=4
mycas1 = mcscf.CASCI(mf,norb,nelec).run()
casci_energy=mycas1.e_tot

# for CASSCF energy
mycas2 = mcscf.CASSCF(mf,norb,nelec).run()
casscf_energy=mycas2.e_tot

# for NEVPT2
ci_nevpt_e1 = mrpt.NEVPT(mycas2).kernel()

with open(file_name,'a') as f:
    np.savetxt(f,[[rhf_energy, casci_energy,casscf_energy,casscf_energy+ci_nevpt_e1]])

