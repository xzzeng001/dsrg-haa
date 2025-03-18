import numpy as np
import openfermion
from openfermion import MolecularData
import scipy
from grouping import ham_grouping 
import os

def prepare_qubit_ham():

    '''
    calculate the effective hamiltonian from the one and two electron integral after the DSRG calculation
    for the following VQE calcualtion
    '''
    os.system('ln -sf ../core_energy.txt .')
    os.system('ln -sf ../oei_a.txt .')
    os.system('ln -sf ../tei_aa.txt .')
    os.system('ln -sf ../tei_ab.txt .')
    
    mm_energy=np.loadtxt("core_energy.txt")
    mm=int(mm_energy[0])
    nuclear_energy=mm_energy[1]
    #print('dims and energy: ',mm,energy)
    
    a=np.loadtxt("oei_a.txt")
    #print(a.reshape(mm,mm,order='C'))
    
    aa=np.loadtxt("tei_aa.txt")
    ab=np.loadtxt("tei_ab.txt")
    #print(ab.reshape(mm,mm,mm,mm,order='C'))
    
    # the geometry is unnecessary
    molecule = MolecularData(geometry="H 0.0 0.0 0.0\nH 0.0 0.0 1.0",
        basis="cc-pvtz",
        multiplicity=1,
        charge=0)
    
    molecule.one_body_integrals = a.reshape(mm,mm,order='C')
    molecule.two_body_integrals = ab.reshape(mm,mm,mm,mm,order='C')-aa.reshape(mm,mm,mm,mm,order='C')
    molecule.nuclear_repulsion  = nuclear_energy
    molecule.save()
    
    ham=molecule.get_molecular_hamiltonian()
    fop=openfermion.transforms.get_fermion_operator(ham)
    qop=openfermion.symmetry_conserving_bravyi_kitaev(fop,8,4)
#    qop = openfermion.jordan_wigner(fop)
    
    np.save('ham.npy',qop)
    
    ham_grouping(qop)
    
    ham_matrix = openfermion.get_sparse_operator(qop)
    e, eigvectors = scipy.sparse.linalg.eigsh(ham_matrix, k=1, which="SA")
    
    with open("fci.txt","a") as f:
       np.savetxt(f,[e])

