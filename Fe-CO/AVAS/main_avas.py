# ===== PySCF AVAS demos for three systems =====
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas

def run_avas(mol, mf, aolabels, ncas_hint=None, nelecas_hint=None, thresh=0.3):
    """
    AVAS wrapper that works across common PySCF APIs.
    Returns: ncore, ncas, nelecas, mo_coeff_av, idx_active (list of MO indices), overlaps (array)
    """
    from pyscf.mcscf import avas
    # API variant 1: avas.AVAS(mf, aolabels, ncas, thresh)
    if ncas_hint is not None and nelecas_hint is not None:
        mycas = mcscf.CASCI(mf, ncas_hint, nelecas_hint)
    else:
        # If not specified, AVAS can suggest ncas from overlap spectrum
        # set a generous upper bound and let threshold prune
        ncas_hint = 8
        nelecas_hint = min(mol.nelectron, 2*ncas_hint)
        mycas = mcscf.CASCI(mf, ncas_hint, nelecas_hint)

    av = avas.AVAS(mycas, aolabels, threshold=thresh)
    ncore, ncas, mo_coeff = av.kernel()
    # Extract active MO indices via overlap spectrum
    ovlp = av.mo_occ     if hasattr(av, 'mo_occ') else None
    idx_active = av.cas_list if hasattr(av, 'cas_list') else list(range(ncore, ncore+ncas))
    return ncore, ncas, mycas.nelecas, mo_coeff, idx_active, ovlp


# ---------- Linear Fe–CO fragment ----------
# Minimal model: Fe–CO along z; keep small basis + ECP for stability.
mol = gto.Mole()
mol.build(atom='feco.xyz', basis={'Fe':'lanl2dz','C':'def2-svp','O':'def2-svp'},
        ecp={'Fe':'lanl2dz'},
        charge=0, spin=0, verbose = 4)

mf = scf.UHF(mol)  # open-shell often helps convergence; switch to RHF if closed-shell
mf.max_cycle = 200
mf.run()

n_elec=sum(mol.nelec)
n_orb=mol.nao_nr()
print('n_elec: ',n_elec)
print('n_orb: ',n_orb)

import time

start_time=time.time()

# Seed AVAS with Fe 3d/4s/4p and CO 2p to capture d–π back-bonding and σ donation
aolabels_fe_co = ['Fe 3d', 'Fe 4s', 'Fe 4p', 'C 2p', 'O 2p']
ncore, ncas, nelecas, moC, idx_act, weights = run_avas(mol, mf, aolabels_fe_co, ncas_hint=8, nelecas_hint=10, thresh=0.2)
end_time=time.time()

print('[Fe–CO] ncore, ncas, nelecas =', ncore, ncas, nelecas)
print('[Fe–CO] active MO idx =', idx_act)
print('Total time (s): ',end_time-start_time)

