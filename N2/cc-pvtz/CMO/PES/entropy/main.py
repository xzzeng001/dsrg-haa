"""Basic example script.

This script is a modified version from the ground state cas calculation
in scine_autocas.main_functions
"""
# -*- coding: utf-8 -*-
__copyright__ = """This file is part of SCINE AutoCAS.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details
"""

import os

from scine_autocas import Autocas
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.main_functions import MainFunctions
from scine_autocas.plots.entanglement_plot import EntanglementPlot

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

file_name='energy_pes.txt'
np.savetxt(file_name,['# bond-length  CASCI  CASSCF  CASSCF+NEVPT2'],fmt='%s')
current_path = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    for rr in np.array([0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0]):
        file_path0 = current_path+'/'+str(rr.round(1))
        file_orb_list=file_path0+'/orb_list.txt'
        os.mkdir(file_path0)

        path_to_this_file = file_path0
        file_xyz=open(path_to_this_file+'/n2.xyz',mode='w')
        file_xyz.writelines(["2\n","n2 \n","N 0.0  0.0  0.0\n","N 0.0  0.0  "+str(rr)])
        file_xyz.close
        file_xyz.flush()

        xyz = path_to_this_file + '/n2.xyz'

        # The main autoCAS process 
        main_functions = MainFunctions()
        new_molecule = Molecule(xyz)
        new_autocas = Autocas(new_molecule)
        new_interface = Molcas([new_molecule])

        new_interface.project_name = "n2"
        new_interface.settings.work_dir = path_to_this_file + "/results/work_dir"
        new_interface.settings.xyz_file = xyz
        new_interface.environment.molcas_scratch_dir = path_to_this_file + "/results/scratch"

        new_interface.settings.method = "dmrg-ci"
        new_interface.settings.basis_set="cc-pvtz"
        new_interface.settings.dmrg_bond_dimension = 300
        new_interface.settings.dmrg_sweeps = 5
        cas_occ, cas_index = main_functions.conventional(new_autocas, new_interface)

        n_elec = sum(cas_occ)
        n_orb = len(cas_occ)

        # for pyscf process
        mol = gto.Mole()
        mol.atom = [["N", [0., 0., 0.0]],["N", [0., 0., rr]]]
        mol.basis = 'cc-pvtz'
        mol.spin = 0
        mol.build()
        mol.verbose = 4

        mf = scf.RHF(mol)
        mf.kernel()
        e_hf=mf.e_tot

        mycas1 = mcscf.CASCI(mf,n_orb,n_elec)
        ## init_cas add 1 for indice
        orb_indice=[i+1 for i in cas_index]

        with open(file_orb_list,'a') as f:
            np.savetxt(f,[orb_indice],fmt='%d')

        mo1 = mycas1.sort_mo(orb_indice)
        #
        e_casci, e_cas, fcivec, mo, mo_energy=mycas1.kernel(mo1)

        mycas2 = mcscf.CASSCF(mf,n_orb,n_elec)
        e_casscf, e_cas, fcivec, mo, mo_energy=mycas2.kernel(mo1)
        #       

        ci_nevpt_e1 = mrpt.NEVPT(mycas2).kernel()    
   
        with open(file_name,'a') as f:
            np.savetxt(f,[[rr, e_casci, e_casscf, e_casscf+ci_nevpt_e1]])  
