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

# ---------- 1) Ethene (C2H4), π/π* AVAS ----------
mol = gto.Mole()
mol.build(atom='c2h4.xyz', basis='def2-svp', charge=0, spin=0,verbose = 4)

mf = scf.RHF(mol).run()

n_elec=sum(mol.nelec)
n_orb=mol.nao_nr()
print('n_elec: ',n_elec)
print('n_orb: ',n_orb)

import time

start_time=time.time()
# Seed AVAS with carbon 2p-like functions for π space
aolabels_ethene = ['C 2p']  # PySCF will match AO labels containing these tokens
ncore, ncas, nelecas, moC, idx_act, weights = run_avas(mol, mf, aolabels_ethene, ncas_hint=2, nelecas_hint=2, thresh=0.2)
end_time=time.time()

print('[Ethene] ncore, ncas, nelecas =', ncore, ncas, nelecas)
print('[Ethene] active MO idx =', idx_act)
print('Total time (s): ',end_time-start_time)

