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
np.savetxt(file_name,['# bond-length   RHF   MP2   CCSD   CISD'],fmt='%s')
current_path = os.getcwd()

for rr in np.array([0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0]):

    # mkdir the work directory
    file_path0 = current_path+'/'+str(rr.round(1))
    file_mo_energy=file_path0+'/mo_energy.txt'
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

    # for the CCSD energy
    mycc=cc.CCSD(mf).run()
    cc_energy=mycc.e_tot

    # for the CISD energy
    myci = ci.CISD(mf).run()
    ci_energy=myci.e_tot
   
    # for mp2 energy
    mymp = mp.MP2(mf).run()
    mp_energy=mymp.e_tot

    with open(file_name,'a') as f:
        np.savetxt(f,[[rr, e_hf, mp_energy, cc_energy, ci_energy]]) 
