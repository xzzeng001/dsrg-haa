# ==> Import Psi4, NumPy, & SciPy <==
import pyscf
from pyscf import gto,scf
import numpy as np
import scipy.linalg as la

def orbital_rotation_ah(mycas,mo_coeff,fcivec,eris,casdm1,casdm2,e_cas):

    mo = mo_coeff
    nmo = mo_coeff.shape[1]
    ncore = mycas.ncore
    ncas = mycas.ncas
    nocc = ncore + ncas
    r0   = None
    ci0  = None
 
    tol=1e-7   
    conv_tol_grad = np.sqrt(tol)
    norm_gorb = norm_gci = -1
    max_cycle_micro = 4
    conv_tol_ddm = conv_tol_grad * 3
    conv = False
    max_stepsize = 0.02
    casdm1_prev = casdm1_last = casdm1

    rota = mycas.rotate_orb_cc(mo, lambda:fcivec, lambda:casdm1, lambda:casdm2,
                                eris, r0, conv_tol_grad*.3, max_stepsize)

    imicro=0
    for u, g_orb, njk, r0 in rota:
        imicro += 1
        norm_gorb = np.linalg.norm(g_orb)
        if imicro == 1:
            norm_gorb0 = norm_gorb
        norm_t = np.linalg.norm(u-np.eye(nmo))

        if imicro >= max_cycle_micro:
            print('micro %2d  |u-1|=%5.3g  |g[o]|=%5.3g',
                      imicro, norm_t, norm_gorb)
            break

        casdm1, casdm2, gci, fcivec = \
                mycas.update_casdm(mo, u, fcivec, e_cas, eris)

        norm_ddm = np.linalg.norm(casdm1 - casdm1_last)
        norm_ddm_micro = np.linalg.norm(casdm1 - casdm1_prev)
        casdm1_prev = casdm1
        if isinstance(gci, np.ndarray):
            norm_gci = np.linalg.norm(gci)
            print('micro %2d  |u-1|=%5.3g  |g[o]|=%5.3g  |g[c]|=%5.3g  |ddm|=%5.3g',
                      imicro, norm_t, norm_gorb, norm_gci, norm_ddm)
        else:
            norm_gci = None
            print('micro %2d  |u-1|=%5.3g  |g[o]|=%5.3g  |g[c]|=%s  |ddm|=%5.3g',
                      imicro, norm_t, norm_gorb, norm_gci, norm_ddm)


        if (norm_t < conv_tol_grad or
            (norm_gorb < conv_tol_grad*.5 and
             (norm_ddm < conv_tol_ddm*.4 or norm_ddm_micro < conv_tol_ddm*.4))):
            break

    rota.close()
    rota = None

    eris = None
    # keep u, g_orb in locals() so that they can be accessed by callback
    u = u.copy()
    g_orb = g_orb.copy()
    mo = mycas.rotate_mo(mo, u)
    eris = mycas.ao2mo(mo)

    return mo
'''
    max_offdiag_u = np.abs(np.triu(u, 1)).max()
    if max_offdiag_u < mycas.small_rot_tol:
        small_rot = True
    else:
        small_rot = False

    from pyscf.mcscf.addons import StateAverageMCSCFSolver

    if not isinstance(mycas, StateAverageMCSCFSolver):
        # The fcivec from builtin FCI solver is a numpy.ndarray
        if not isinstance(fcivec, np.ndarray):
            fcivec = small_rot
    else:
        newvecs = []
        for subvec in fcivec:
            # CI vector obtained by builtin FCI is a numpy array
            if not isinstance(subvec, np.ndarray):
                newvecs.append(small_rot)
            else:
                newvecs.append(subvec)
        fcivec = newvecs

    return mo, fcivec, eris
#    e_tot, e_cas, fcivec = mycas.casci(mo, fcivec, eris)
'''
