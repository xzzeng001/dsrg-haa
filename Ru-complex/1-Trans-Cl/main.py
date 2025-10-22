import numpy as np
import pyscf
from pyscf import gto,lo,mp
from pyscf import scf,cc,df
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

XYZ_FILE = "Trans-Cl.xyz"
# ===== 1) 分子与 SCF =====
mol = gto.Mole()
mol.build(
    atom=XYZ_FILE,   # 从 xyz 读取
    unit="Angstrom",
    basis="def2-tzvp",
    ecp={"Ru": "ECP28MDF"},            # Ru 用 Stuttgart/Cologne ECP28MDF
    charge=0,
    spin=0,                            # 2S，单重态=0
    symmetry=True,
    symmetry_subgroup="C2v",
    verbose=4
)

mf = scf.RHF(mol)
mf.max_cycle = 200
mf.conv_tol  = 1e-10
mf.kernel()
e_hf = mf.e_tot

mf.mol.incore_anyway = True

nelec = mol.nelectron
norb  = mol.nao_nr()
print(f"Total electrons: {nelec}")
print(f"Total AO orbitals: {norb}")

# 保存 MO 能量
np.savetxt('mo_energy.txt', mf.mo_energy.reshape(1, -1))

# ===== 2) 选活性轨道（示例：HOMO-LUMO 两轨道）=====
occ = mf.mo_occ
homo = np.where(occ > 1e-8)[0].max()   # 最后一个占据轨道索引
lumo = homo + 1
caslst = [homo, lumo]                  # 你也可以改成自己的列表
ncas = len(caslst)
nelecas = int(round(occ[caslst].sum()))  # 活性电子数（RHF 情况下为占据数之和）

print("HOMO index:", homo, "LUMO index:", lumo)
print("CAS list (0-based):", caslst, "nelecas:", nelecas)

# ===== 3) 基于阈值自动冻结轨道（绝对能量大于阈值）=====
E_THRESH = 1.0  # Hartree；可按需改，例如 0.5 或 1.5 Hartree

# 仅用于查看：对能量排序（不改变 MO 排序）
order = np.argsort(mf.mo_energy)
#print("Lowest 10 MO energies (Eh):", mf.mo_energy[order][:10])

# 生成冻结列表：|epsilon| > E_THRESH 的全部轨道
frozen_raw = [i for i, e in enumerate(mf.mo_energy) if abs(e) > E_THRESH]

# 不要冻结活性轨道；并做去重、排序
frozen_set = set(frozen_raw) - set(caslst)
frozen_lst = sorted(frozen_set)

print(f"Energy threshold = {E_THRESH:.3f} Eh")
print(f"#frozen = {len(frozen_lst)}")
np.savetxt("frozen_indices.txt", np.array(frozen_lst, dtype=int), fmt="%d")

# 额外：给出未冻结的虚轨道索引（便于你后续检查/使用）
all_idx   = list(range(mf.mo_coeff.shape[1]))
vir_index = [i for i in all_idx
             if (i not in frozen_lst) and (i not in caslst) and (occ[i] < 1e-6)]
np.savetxt("vir_indices_unfrozen.txt", np.array(vir_index, dtype=int), fmt="%d")

#mo = mf.mo_coeff.copy()
# 把 MO 排成 [核心 | 活性] 的顺序（非常重要！CASCI按这个切片）
#mo_reduced = np.hstack((mo[:, :ncore], mo[:, caslst]))  # shape: (nao, 95)

# ===== 4) 运行 CASSCF（冻结 frozen_lst）=====
mc = mcscf.CASSCF(mf, ncas, nelecas)
mc.frozen = frozen_lst               # 关键：冻结阈值外轨道
# 把 caslst 放入活性空间（索引按当前 MO 排序）
mo_init = mc.sort_mo(caslst, base=0)
mc.max_cycle_macro = 1
mc.conv_tol = 1e-9
mc.max_memory  = 64000
mc.fcisolver.max_memory = 64000

e_tot, e_cas, fcivec, mo_coeff, mo_energy=mc.kernel(mo_init)

base_e_corr=e_tot-e_hf

print("HF energy (Eh):", e_hf)
print("CASSCF energy (Eh):", mc.e_tot)

mbe_fci_corr(mf,e_tot,base_e_corr,homo,ncas,nelecas,caslst,vir_index,frozen_lst,mo_coeff,2)

import sys
sys.exit(0)

all_active_lst=[]
ind1,ee,ind1_corr=np.loadtxt('one-body.txt',delimiter=",",unpack=True)
ind2,ind3,ee_2,ind2_corr=np.loadtxt('two-body.txt',delimiter=",",unpack=True)
tol=0.002

all_active_lst=copy.copy(caslst)
for ii,rr in enumerate(ind1_corr):
    if abs(rr) > tol and int(ind1[ii]) not in all_active_lst:
        all_active_lst.append(int(ind1[ii]))

for ii,rr in enumerate(ind2_corr):
    if abs(rr) > tol:
        if int(ind2[ii]) not in all_active_lst:
            all_active_lst.append(int(ind2[ii]))
        if int(ind3[ii]) not in caslst:
            all_active_lst.append(int(ind3[ii]))

n_all_active_orb=len(all_active_lst)

n_all_active_elec=0
for ii in all_active_lst:
    if ii < ind_homo:
        n_all_active_elec += 2

#######################################
# generate hamilonian
#######################################
mycas1 = mcscf.CASCI(mf,n_all_active_orb,n_all_active_elec)
mycas1.max_cycle = 500
mycas1.frozen=frozen_lst

orb_indice=[i+1 for i in all_active_lst]
mo1 = mycas1.sort_mo(orb_indice,mo_coeff=mo_coeff)
e_tot, e_cas, fcivec, mo_coeff, mo_energy=mycas1.kernel(mo1)

generate_ham(mycas1,mo_coeff,mycas1.ncas)

print('energy CASCI: ',e_tot)
base_e_corr=e_tot-e_hf
print('base correlation energy :', base_e_corr)
sys.stdout.flush()


'''
#######################################
# orbital roation
#######################################
ci0=None
eris = mycas.ao2mo(mo_coeff)      
e_tot, e_cas, fcivec = mycas.casci(mo_coeff, ci0, eris)
casdm1, casdm2 = mycas.fcisolver.make_rdm12(fcivec, mycas.ncas, mycas.nelecas)
mo_coeff=orbital_rotation_ah(mycas,mf.mo_coeff,fcivec,eris,casdm1,casdm2,e_cas)
e_tot, e_cas, fcivec=mycas.casci(mo_coeff)

generate_ham(mycas,mo_coeff,mycas.ncas)

print('energy CASCI: ',e_tot)
base_e_corr=e_tot-e_hf

print('base correlation energy :', base_e_corr)
sys.stdout.flush()

########################################
# The NEVPT2 for left orbitals
########################################
ci_nevpt_e1 = mrpt.NEVPT(mycas).kernel()
print('nevpt2: ',ci_nevpt_e1)
sys.stdout.flush()
'''
