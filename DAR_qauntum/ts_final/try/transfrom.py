import numpy as np
import openfermion
from openfermion import MolecularData
import scipy

mm_energy=np.loadtxt("core_energy.txt")
mm=int(mm_energy[0])
energy=mm_energy[1]
#print('dims and energy: ',mm,energy)

a=np.loadtxt("oei_a.txt")
#print(a.reshape(mm,mm,order='C'))

aa=np.loadtxt("tei_aa.txt")
ab=np.loadtxt("tei_ab.txt")
#print(ab.reshape(mm,mm,mm,mm,order='C'))

molecule = MolecularData(geometry="H 0.0 0.0 0.0\nF 0.0 0.0 1.5",
    basis="cc-pvdz",
    multiplicity=1,
    charge=0)

molecule.one_body_integrals = a.reshape(mm,mm,order='C')
molecule.two_body_integrals = ab.reshape(mm,mm,mm,mm,order='C')-aa.reshape(mm,mm,mm,mm,order='C')
molecule.nuclear_repulsion  = energy
molecule.save()

ham=molecule.get_molecular_hamiltonian()
fop=openfermion.transforms.get_fermion_operator(ham)
qop = openfermion.jordan_wigner(fop)

ham_matrix = openfermion.get_sparse_operator(qop)
e, eigvectors = scipy.sparse.linalg.eigsh(ham_matrix, k=1, which="SA")

with open("energy.txt","a") as f:
    np.savetxt(f,[e])

print("e: ",e)

###print("H: ",qop)
