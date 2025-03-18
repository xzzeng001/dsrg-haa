import numpy as np
import openfermion
from openfermion import MolecularData

def generate_ham(mycas,mo,n_orb):
    ########################################
    # generate effective hamiltonian
    ########################################

    one_body_mo,energy_core=mycas.get_h1eff(mo)
    two_electron_compressed=mycas.get_h2eff(mo)
    two_electron_integrals = np.empty((n_orb, n_orb,n_orb, n_orb))
    
    # Unpack symmetry.
    n_pairs = n_orb * (n_orb + 1) // 2
    if two_electron_compressed.ndim == 2:
    
        # Case of 4-fold symmetry.
        assert(two_electron_compressed.size == n_pairs ** 2)
        pq = 0
        for p in range(n_orb):
            for q in range(p + 1):
                rs = 0
                for r in range(n_orb):
                    for s in range(r + 1):
                        pqrs_value = two_electron_compressed[pq, rs]
                        two_electron_integrals[p, s, r, q] = float(pqrs_value)
                        two_electron_integrals[q, s, r, p] = float(pqrs_value)
                        two_electron_integrals[p, r, s, q] = float(pqrs_value)
                        two_electron_integrals[q, r, s, p] = float(pqrs_value)
                        rs += 1
                pq += 1
    else:
    
        # Case of 8-fold symmetry.
        assert(two_electron_compressed.size == n_pairs * (n_pairs + 1) // 2)
        pq = 0
        pqrs = 0
        for p in range(n_orb):
            for q in range(p + 1):
                rs = 0
                for r in range(p + 1):
                    for s in range(r + 1):
                        if pq >= rs:
                            pqrs_value = two_electron_compressed[pqrs]
                            two_electron_integrals[p, s,
                                                   r, q] = float(pqrs_value)
                            two_electron_integrals[q, s,
                                                   r, p] = float(pqrs_value)
                            two_electron_integrals[p, r,
                                                   s, q] = float(pqrs_value)
                            two_electron_integrals[q, r,
                                                   s, p] = float(pqrs_value)
                            two_electron_integrals[s, p,
                                                   q, r] = float(pqrs_value)
                            two_electron_integrals[s, q,
                                                   p, r] = float(pqrs_value)
                            two_electron_integrals[r, p,
                                                   q, s] = float(pqrs_value)
                            two_electron_integrals[r, q,
                                                   p, s] = float(pqrs_value)
                            pqrs += 1
                        rs += 1
                pq += 1
    
    two_body_mo=two_electron_integrals
    
    #
    molecule = MolecularData(geometry="H 0.0 0.0 0.0\nH 0.0 0.0 0.7",
        basis="sto-3g",
        multiplicity=1,
        charge=0)
    
    molecule.one_body_integrals = one_body_mo
    molecule.two_body_integrals = two_body_mo
    molecule.nuclear_repulsion  = energy_core
    #molecule.save()
    
    ham=molecule.get_molecular_hamiltonian()
    hamiltonian_fermOp=openfermion.transforms.get_fermion_operator(ham)
    
    
    n_qubits = openfermion.count_qubits(hamiltonian_fermOp)
    
    print('all qubits: ',n_qubits)
    
    hamiltonian_qubitOp = openfermion.jordan_wigner(hamiltonian_fermOp)
    
    np.save('ham.npy',hamiltonian_qubitOp)
    print('ham saved!!!')
    
    
