import numpy as np
from pyscf import mcscf,mrpt
import copy
import sys

def mbe_fci_corr(mf,e_tot,base_e_corr,ind_homo,ncas,nelec,caslst,vir_index,frozen_lst,mo_coeff,norder,file_one_body,file_two_body):

    ## for the one-orbital added
    if norder == 1:

        nl_occ=len(vir_index)
        one_e_corr=np.zeros(nl_occ)
        one_e_tot=np.zeros(nl_occ)
        for ii, kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
            new_cas=copy.copy(caslst)
            new_cas.append(kk)
            ncas=len(new_cas)
            mycas1 = mcscf.CASCI(mf,ncas,new_nelec1)
            mycas1.max_cycle = 500
            mycas1.frozen=frozen_lst
            orb_indice1=[i+1 for i in new_cas]
            mo1 = mycas1.sort_mo(orb_indice1,mo_coeff=mo_coeff)
            e_tot1=mycas1.kernel(mo1)[0]
            one_e_tot[ii]=e_tot1
            one_e_corr[ii]=e_tot1-e_tot
            with open(file_one_body,'a') as f:
                np.savetxt(f,[[kk,e_tot1,e_tot1-e_tot]],fmt='%d,%f,%f')
                
            del new_cas, mycas1
        
        print('one-all :', base_e_corr+np.sum(one_e_corr))
        sys.stdout.flush()
    
    elif norder == 2:

        # for the one-body correlation
        nl_occ=len(vir_index)
        one_e_corr=np.zeros(nl_occ)
        one_e_tot=np.zeros(nl_occ)
        for ii, kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
            new_cas=copy.copy(caslst)
            new_cas.append(kk)
            ncas=len(new_cas)
            mycas1 = mcscf.CASCI(mf,ncas,new_nelec1)
            mycas1.max_cycle = 500
            mycas1.frozen=frozen_lst
            orb_indice1=[i+1 for i in new_cas]
            mo1 = mycas1.sort_mo(orb_indice1,mo_coeff=mo_coeff)
            e_tot1=mycas1.kernel(mo1)[0]
            one_e_tot[ii]=e_tot1
            one_e_corr[ii]=e_tot1-e_tot
            with open(file_one_body,'a') as f:
                np.savetxt(f,[[kk,e_tot1,e_tot1-e_tot]],fmt='%d,%f,%f')
                
            del new_cas, mycas1
        
        print('one-all :', base_e_corr+np.sum(one_e_corr))
        sys.stdout.flush()
        
        ## for the two-orbital added
        two_e_corr=np.zeros((nl_occ,nl_occ))
        two_e_tot=np.zeros((nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
           
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    new_cas=copy.copy(caslst)
                    new_cas.append(kk)
                    new_cas.append(ll)
                    ncas=len(new_cas)
                    mycas2 = mcscf.CASCI(mf,ncas,new_nelec2)
                    orb_indice1=[i+1 for i in new_cas]
                    mo1 = mycas2.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                    mycas2.max_cycle = 500
                    mycas2.frozen=frozen_lst
                    e_tot1=mycas2.kernel(mo1)[0]
                    two_e_tot[ii][jj]=e_tot1
                    two_e_corr[ii][jj]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]
                    with open(file_two_body,'a') as f:
                        np.savetxt(f,[[kk,ll,e_tot1,two_e_corr[ii][jj]]],fmt='%d,%d,%f,%f')
        
                    del new_cas, mycas2
        
        print('two-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr))
        sys.stdout.flush()
    
    else:
 
        # for the one-body correlation
        nl_occ=len(vir_index)
        one_e_corr=np.zeros(nl_occ)
        one_e_tot=np.zeros(nl_occ)
        for ii, kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
            new_cas=copy.copy(caslst)
            new_cas.append(kk)
            ncas=len(new_cas)
            mycas1 = mcscf.CASCI(mf,ncas,new_nelec1)
            mycas1.max_cycle = 500
            mycas1.frozen=frozen_lst
            orb_indice1=[i+1 for i in new_cas]
            mo1 = mycas1.sort_mo(orb_indice1,mo_coeff=mo_coeff)
            e_tot1=mycas1.kernel(mo1)[0]
            one_e_tot[ii]=e_tot1
            one_e_corr[ii]=e_tot1-e_tot
            with open('one-body.txt','a') as f:
                np.savetxt(f,[[kk,e_tot1,e_tot1-e_tot]],fmt='%d,%f,%f')
                
            del new_cas, mycas1
        
        print('one-all :', base_e_corr+np.sum(one_e_corr))
        sys.stdout.flush()
        
        ## for the two-orbital added
        two_e_corr=np.zeros((nl_occ,nl_occ))
        two_e_tot=np.zeros((nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
           
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    new_cas=copy.copy(caslst)
                    new_cas.append(kk)
                    new_cas.append(ll)
                    ncas=len(new_cas)
                    mycas2 = mcscf.CASCI(mf,ncas,new_nelec2)
                    orb_indice1=[i+1 for i in new_cas]
                    mo1 = mycas2.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                    mycas2.max_cycle = 500
                    mycas2.frozen=frozen_lst
                    e_tot1=mycas2.kernel(mo1)[0]
                    two_e_tot[ii][jj]=e_tot1
                    two_e_corr[ii][jj]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]
                    with open('two-body.txt','a') as f:
                        np.savetxt(f,[[kk,ll,e_tot1,two_e_corr[ii][jj]]],fmt='%d,%d,%f,%f')
        
                    del new_cas, mycas2
        
        print('two-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr))
        sys.stdout.flush()
   
        ## for the three-orbital added
        three_e_corr=np.zeros((nl_occ,nl_occ,nl_occ))
        three_e_tot=np.zeros((nl_occ,nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
        
            for jj,ll in enumerate(vir_index):
                if ll > kk:            
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    for ij,lk in enumerate(vir_index):
                        if lk > ll:
                            if lk < ind_homo:
                                new_nelec3=new_nelec2+2
                            else:
                                new_nelec3=new_nelec2
        
                            new_cas=copy.copy(caslst)
                            new_cas.append(kk)
                            new_cas.append(ll)
                            new_cas.append(lk)
                            ncas=len(new_cas)
                            mycas3 = mcscf.CASCI(mf,ncas,new_nelec3)
                            mycas3.max_cycle = 500
                            mycas3.frozen=frozen_lst
        
                            orb_indice1=[i+1 for i in new_cas]
                            mo1 = mycas3.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                            e_tot1=mycas3.kernel(mo1)[0]
                            three_e_tot[ii][jj][ij]=e_tot1
                            three_e_corr[ii][jj][ij]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]-one_e_corr[ij]-two_e_corr[ii][jj]-two_e_corr[ii][ij]-two_e_corr[jj][ij]
                            with open('three-body.txt','a') as f:
                                np.savetxt(f,[[kk,ll,lk,e_tot1,three_e_corr[ii][jj][ij]]],fmt='%d,%d,%d,%f,%f')
        
                            del new_cas, mycas3
        
        print('three-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr)+np.sum(three_e_corr))
        sys.stdout.flush()
        
        ## for the four-orbital added
        four_e_corr=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ))
        four_e_tot=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
        
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    for ij,lk in enumerate(vir_index):
                        if lk > ll: 
                            if lk < ind_homo:
                                new_nelec3=new_nelec2+2
                            else:
                                new_nelec3=new_nelec2
        
                            for im,kn in enumerate(vir_index):
                                if  kn > lk:
                                    if kn < ind_homo:
                                        new_nelec4=new_nelec3+2
                                    else:
                                        new_nelec4=new_nelec3
        
                                    new_cas=copy.copy(caslst)
                                    new_cas.append(kk)
                                    new_cas.append(ll)
                                    new_cas.append(lk)
                                    new_cas.append(kn)
                                    ncas=len(new_cas)
                                    mycas4 = mcscf.CASCI(mf,ncas,new_nelec4)
                                    mycas4.max_cycle = 500
                                    mycas4.frozen=frozen_lst
        
                                    orb_indice1=[i+1 for i in new_cas]
                                    mo1 = mycas4.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                                    e_tot1=mycas4.kernel(mo1)[0]
                                    four_e_tot[ii][jj][ij][im]=e_tot1
                                    four_e_corr[ii][jj][ij][im]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]-one_e_corr[ij]-one_e_corr[im]-two_e_corr[ii][jj]-two_e_corr[ii][ij]-two_e_corr[ii][im]-two_e_corr[jj][ij]-two_e_corr[jj][im]-two_e_corr[ij][im]-three_e_corr[ii][jj][ij]-three_e_corr[ii][jj][im]-three_e_corr[ii][ij][im]-three_e_corr[jj][ij][im]
                                    with open('four-body.txt','a') as f:
                                        np.savetxt(f,[[kk,ll,lk,kn,e_tot1,four_e_corr[ii][jj][ij][im]]],fmt='%d,%d,%d,%d,%f,%f')
        
                                    del new_cas, mycas4
        
        print('four-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr)+np.sum(three_e_corr)+np.sum(four_e_corr))
        sys.stdout.flush()
        
        ## for the five-orbital added
        five_e_corr=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        five_e_tot=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
        
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    for ij,lk in enumerate(vir_index):
                        if lk > ll:
                            if lk < ind_homo:
                                new_nelec3=new_nelec2+2
                            else:
                                new_nelec3=new_nelec2
        
                            for im,kn in enumerate(vir_index):
                                if kn > lk:
                                    if kn < ind_homo:
                                        new_nelec4=new_nelec3+2
                                    else:
                                        new_nelec4=new_nelec3
        
                                    for io,ko in enumerate(vir_index):
                                        if ko > kn:
                                            if ko < ind_homo:
                                                new_nelec5=new_nelec4+2
                                            else:
                                                new_nelec5=new_nelec4
        
                                            new_cas=copy.copy(caslst)
                                            new_cas.append(kk)
                                            new_cas.append(ll)
                                            new_cas.append(lk)
                                            new_cas.append(kn)
                                            new_cas.append(ko)
                                            ncas=len(new_cas)
                                            mycas5 = mcscf.CASCI(mf,ncas,new_nelec5)
                                            mycas5.max_cycle = 500
                                            mycas5.frozen=frozen_lst
        
                                            orb_indice1=[i+1 for i in new_cas]
                                            mo1 = mycas5.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                                            e_tot1=mycas5.kernel(mo1)[0]
                                            five_e_tot[ii][jj][ij][im][io]=e_tot1
                                            five_e_corr[ii][jj][ij][im][io]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]-one_e_corr[ij]-one_e_corr[im]-one_e_corr[io]-two_e_corr[ii][jj]-two_e_corr[ii][ij]-two_e_corr[ii][im]-two_e_corr[ii][io]-two_e_corr[jj][ij]-two_e_corr[jj][im]-two_e_corr[jj][io]-two_e_corr[ij][im]-two_e_corr[ij][io]-two_e_corr[im][io]-three_e_corr[ii][jj][ij]-three_e_corr[ii][jj][im]-three_e_corr[ii][jj][io]-three_e_corr[ii][ij][im]-three_e_corr[ii][ij][io]-three_e_corr[ii][im][io]-three_e_corr[jj][ij][im]-three_e_corr[jj][ij][io]-three_e_corr[jj][im][io]-three_e_corr[ij][im][io]-four_e_corr[ii][jj][ij][im]-four_e_corr[ii][jj][ij][io]-four_e_corr[ii][jj][im][io]-four_e_corr[ii][ij][im][io]-four_e_corr[jj][ij][im][io]
                                            with open('five-body.txt','a') as f:
                                                np.savetxt(f,[[kk,ll,lk,kn,ko,e_tot1,five_e_corr[ii][jj][ij][im][io]]],fmt='%d,%d,%d,%d,%d,%f,%f')
        
                                            del new_cas, mycas5
        
        print('five-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr)+np.sum(three_e_corr)+np.sum(four_e_corr)+np.sum(five_e_corr))
        sys.stdout.flush()
        
        ## for the six-orbital added
        six_e_corr=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        six_e_tot=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
        
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    for ij,lk in enumerate(vir_index):
                        if lk > ll:
                            if lk < ind_homo:
                                new_nelec3=new_nelec2+2
                            else:
                                new_nelec3=new_nelec2
        
                            for im,kn in enumerate(vir_index):
                                if kn > lk:
                                    if kn < ind_homo:
                                        new_nelec4=new_nelec3+2
                                    else:
                                        new_nelec4=new_nelec3
        
                                    for io,ko in enumerate(vir_index):
                                        if ko > kn:
                                            if ko < ind_homo:
                                                new_nelec5=new_nelec4+2
                                            else:
                                                new_nelec5=new_nelec4
        
                                            for i6,k6 in enumerate(vir_index):
                                                if k6 > ko:
                                                    if k6 < ind_homo:
                                                        new_nelec6=new_nelec5+2
                                                    else:
                                                        new_nelec6=new_nelec5
        
        
                                                    new_cas=copy.copy(caslst)
                                                    new_cas.append(kk)
                                                    new_cas.append(ll)
                                                    new_cas.append(lk)
                                                    new_cas.append(kn)
                                                    new_cas.append(ko)
                                                    new_cas.append(k6)
                                                    ncas=len(new_cas)
                                                    mycas6 = mcscf.CASCI(mf,ncas,new_nelec6)
                                                    mycas6.max_cycle = 500
                                                    mycas6.frozen=frozen_lst
        
                                                    orb_indice1=[i+1 for i in new_cas]
                                                    mo1 = mycas6.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                                                    e_tot1=mycas6.kernel(mo1)[0]
                                                    six_e_tot[ii][jj][ij][im][io][i6]=e_tot1
                                                    six_e_corr[ii][jj][ij][im][io][i6]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]-one_e_corr[ij]-one_e_corr[im]-one_e_corr[io]-one_e_corr[i6]-two_e_corr[ii][jj]-two_e_corr[ii][ij]-two_e_corr[ii][im]-two_e_corr[ii][io]-two_e_corr[ii][i6]-two_e_corr[jj][ij]-two_e_corr[jj][im]-two_e_corr[jj][io]-two_e_corr[jj][i6]-two_e_corr[ij][im]-two_e_corr[ij][io]-two_e_corr[ij][i6]-two_e_corr[im][io]-two_e_corr[im][i6]-two_e_corr[io][i6]-three_e_corr[ii][jj][ij]-three_e_corr[ii][jj][im]-three_e_corr[ii][jj][io]-three_e_corr[ii][jj][i6]-three_e_corr[ii][ij][im]-three_e_corr[ii][ij][io]-three_e_corr[ii][ij][i6]-three_e_corr[ii][im][io]-three_e_corr[ii][im][i6]-three_e_corr[ii][io][i6]-three_e_corr[jj][ij][im]-three_e_corr[jj][ij][io]-three_e_corr[jj][ij][i6]-three_e_corr[jj][im][io]-three_e_corr[jj][im][i6]-three_e_corr[jj][io][i6]-three_e_corr[ij][im][io]-three_e_corr[ij][im][i6]-three_e_corr[ij][io][i6]-three_e_corr[im][io][i6]-four_e_corr[ii][jj][ij][im]-four_e_corr[ii][jj][ij][io]-four_e_corr[ii][jj][ij][i6]-four_e_corr[ii][jj][im][io]-four_e_corr[ii][jj][im][i6]-four_e_corr[ii][jj][io][i6]-four_e_corr[ii][ij][im][io]-four_e_corr[ii][ij][im][i6]-four_e_corr[ii][ij][io][i6]-four_e_corr[ii][im][io][i6]-four_e_corr[jj][ij][im][io]-four_e_corr[jj][ij][im][i6]-four_e_corr[jj][ij][io][i6]-four_e_corr[jj][im][io][i6]-four_e_corr[ij][im][io][i6]-five_e_corr[ii][jj][ij][im][io]-five_e_corr[ii][jj][ij][im][i6]-five_e_corr[ii][jj][ij][io][i6]-five_e_corr[ii][jj][im][io][i6]-five_e_corr[ii][ij][im][io][i6]-five_e_corr[jj][ij][im][io][i6]
                                                    with open('six-body.txt','a') as f:
                                                        np.savetxt(f,[[kk,ll,lk,kn,ko,k6,e_tot1,six_e_corr[ii][jj][ij][im][io][i6]]],fmt='%d,%d,%d,%d,%d,%d,%f,%f')
                                                    del new_cas, mycas6
        
        print('six-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr)+np.sum(three_e_corr)+np.sum(four_e_corr)+np.sum(five_e_corr)+np.sum(six_e_corr))
        sys.stdout.flush()
        
        ## for the seven-orbital added
        seven_e_corr=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        seven_e_tot=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
        
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    for ij,lk in enumerate(vir_index):
                        if lk > ll:
                            if lk < ind_homo:
                                new_nelec3=new_nelec2+2
                            else:
                                new_nelec3=new_nelec2
        
                            for im,kn in enumerate(vir_index):
                                if kn > lk:
                                    if kn < ind_homo:
                                        new_nelec4=new_nelec3+2
                                    else:
                                        new_nelec4=new_nelec3
        
                                    for io,ko in enumerate(vir_index):
                                        if ko > kn:
                                            if ko < ind_homo:
                                                new_nelec5=new_nelec4+2
                                            else:
                                                new_nelec5=new_nelec4
                                                
                                            for i6,k6 in enumerate(vir_index):
                                                if k6 > ko:
                                                    if k6 < ind_homo:
                                                        new_nelec6=new_nelec5+2
                                                    else:
                                                        new_nelec6=new_nelec5
        
                                                    for i7,k7 in enumerate(vir_index):
                                                        if k7 > k6:
                                                            if k7 < ind_homo:
                                                                new_nelec7=new_nelec6+2
                                                            else:
                                                                new_nelec7=new_nelec6
        
        
                                                            new_cas=copy.copy(caslst)
                                                            new_cas.append(kk)
                                                            new_cas.append(ll)
                                                            new_cas.append(lk)
                                                            new_cas.append(kn)
                                                            new_cas.append(ko)
                                                            new_cas.append(k6)
                                                            new_cas.append(k7)
                                                            ncas=len(new_cas)
                                                            mycas7 = mcscf.CASCI(mf,ncas,new_nelec7)
                                                            mycas7.max_cycle = 500
                                                            mycas7.frozen=frozen_lst
        
                                                            orb_indice1=[i+1 for i in new_cas]
                                                            mo1 = mycas7.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                                                            e_tot1=mycas7.kernel(mo1)[0]
                                                            seven_e_tot[ii][jj][ij][im][io][i6][i7]=e_tot1
                                                            seven_e_corr[ii][jj][ij][im][io][i6][i7]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]-one_e_corr[ij]-one_e_corr[im]-one_e_corr[io]-one_e_corr[i6]-one_e_corr[i7]-two_e_corr[ii][jj]-two_e_corr[ii][ij]-two_e_corr[ii][im]-two_e_corr[ii][io]-two_e_corr[ii][i6]-two_e_corr[ii][i7]-two_e_corr[jj][ij]-two_e_corr[jj][im]-two_e_corr[jj][io]-two_e_corr[jj][i6]-two_e_corr[jj][i7]-two_e_corr[ij][im]-two_e_corr[ij][io]-two_e_corr[ij][i6]-two_e_corr[ij][i7]-two_e_corr[im][io]-two_e_corr[im][i6]-two_e_corr[im][i7]-two_e_corr[io][i6]-two_e_corr[io][i7]-two_e_corr[i6][i7]-three_e_corr[ii][jj][ij]-three_e_corr[ii][jj][im]-three_e_corr[ii][jj][io]-three_e_corr[ii][jj][i6]-three_e_corr[ii][jj][i7]-three_e_corr[ii][ij][im]-three_e_corr[ii][ij][io]-three_e_corr[ii][ij][i6]-three_e_corr[ii][ij][i7]-three_e_corr[ii][im][io]-three_e_corr[ii][im][i6]-three_e_corr[ii][im][i7]-three_e_corr[ii][io][i6]-three_e_corr[ii][io][i7]-three_e_corr[ii][i6][i7]-three_e_corr[jj][ij][im]-three_e_corr[jj][ij][io]-three_e_corr[jj][ij][i6]-three_e_corr[jj][ij][i7]-three_e_corr[jj][im][io]-three_e_corr[jj][im][i6]-three_e_corr[jj][im][i7]-three_e_corr[jj][io][i6]-three_e_corr[jj][io][i7]-three_e_corr[jj][i6][i7]-three_e_corr[ij][im][io]-three_e_corr[ij][im][i6]-three_e_corr[ij][im][i7]-three_e_corr[ij][io][i6]-three_e_corr[ij][io][i7]-three_e_corr[ij][i6][i7]-three_e_corr[im][io][i6]-three_e_corr[im][io][i7]-three_e_corr[im][i6][i7]-three_e_corr[io][i6][i7]-four_e_corr[ii][jj][ij][im]-four_e_corr[ii][jj][ij][io]-four_e_corr[ii][jj][ij][i6]-four_e_corr[ii][jj][ij][i7]-four_e_corr[ii][jj][im][io]-four_e_corr[ii][jj][im][i6]-four_e_corr[ii][jj][im][i7]-four_e_corr[ii][jj][io][i6]-four_e_corr[ii][jj][io][i7]-four_e_corr[ii][jj][i6][i7]-four_e_corr[ii][ij][im][io]-four_e_corr[ii][ij][im][i6]-four_e_corr[ii][ij][im][i7]-four_e_corr[ii][ij][io][i6]-four_e_corr[ii][ij][io][i7]-four_e_corr[ii][ij][i6][i7]-four_e_corr[ii][im][io][i6]-four_e_corr[ii][im][io][i7]-four_e_corr[ii][im][i6][i7]-four_e_corr[ii][io][i6][i7]-four_e_corr[jj][ij][im][io]-four_e_corr[jj][ij][im][i6]-four_e_corr[jj][ij][im][i7]-four_e_corr[jj][ij][io][i6]-four_e_corr[jj][ij][io][i7]-four_e_corr[jj][ij][i6][i7]-four_e_corr[jj][im][io][i6]-four_e_corr[jj][im][io][i7]-four_e_corr[jj][im][i6][i7]-four_e_corr[jj][io][i6][i7]-four_e_corr[ij][im][io][i6]-four_e_corr[ij][im][io][i7]-four_e_corr[ij][im][i6][i7]-four_e_corr[ij][io][i6][i7]-four_e_corr[im][io][i6][i7]-five_e_corr[ii][jj][ij][im][io]-five_e_corr[ii][jj][ij][im][i6]-five_e_corr[ii][jj][ij][im][i7]-five_e_corr[ii][jj][ij][io][i6]-five_e_corr[ii][jj][ij][io][i7]-five_e_corr[ii][jj][ij][i6][i7]-five_e_corr[ii][jj][im][io][i6]-five_e_corr[ii][jj][im][io][i7]-five_e_corr[ii][jj][im][i6][i7]-five_e_corr[ii][jj][io][i6][i7]-five_e_corr[ii][ij][im][io][i6]-five_e_corr[ii][ij][im][io][i7]-five_e_corr[ii][ij][im][i6][i7]-five_e_corr[ii][ij][io][i6][i7]-five_e_corr[ii][im][io][i6][i7]-five_e_corr[jj][ij][im][io][i6]-five_e_corr[jj][ij][im][io][i7]-five_e_corr[jj][ij][im][i6][i7]-five_e_corr[jj][ij][io][i6][i7]-five_e_corr[jj][im][io][i6][i7]-five_e_corr[ij][im][io][i6][i7]-six_e_corr[ii][jj][ij][im][io][i6]-six_e_corr[ii][jj][ij][im][io][i7]-six_e_corr[ii][jj][ij][im][i6][i7]-six_e_corr[ii][jj][ij][io][i6][i7]-six_e_corr[ii][jj][im][io][i6][i7]-six_e_corr[ii][ij][im][io][i6][i7]-six_e_corr[jj][ij][im][io][i6][i7]
                                                            with open('seven-body.txt','a') as f:
                                                                np.savetxt(f,[[kk,ll,lk,kn,ko,k6,k7,e_tot1,seven_e_corr[ii][jj][ij][im][io][i6][i7]]],fmt='%d,%d,%d,%d,%d,%d,%d,%f,%f')
                                                            del new_cas, mycas7
        
        print('seven-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr)+np.sum(three_e_corr)+np.sum(four_e_corr)+np.sum(five_e_corr)+np.sum(six_e_corr)+np.sum(seven_e_corr))
        sys.stdout.flush()
        
        ## for the eight-orbital added
        eight_e_corr=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        eight_e_tot=np.zeros((nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ,nl_occ))
        for ii,kk in enumerate(vir_index):
            if kk < ind_homo:
                new_nelec1=nelec+2
            else:
                new_nelec1=nelec
        
            for jj,ll in enumerate(vir_index):
                if ll > kk:
                    if ll < ind_homo:
                        new_nelec2=new_nelec1+2
                    else:
                        new_nelec2=new_nelec1
        
                    for ij,lk in enumerate(vir_index):
                        if lk > ll:
                            if lk < ind_homo:
                                new_nelec3=new_nelec2+2
                            else:
                                new_nelec3=new_nelec2
        
                            for im,kn in enumerate(vir_index):
                                if kn > lk:
                                    if kn < ind_homo:
                                        new_nelec4=new_nelec3+2
                                    else:
                                        new_nelec4=new_nelec3
        
                                    for io,ko in enumerate(vir_index):
                                        if ko > kn:
                                            if ko < ind_homo:
                                                new_nelec5=new_nelec4+2
                                            else:
                                                new_nelec5=new_nelec4
        
                                            for i6,k6 in enumerate(vir_index):
                                                if k6 > ko:
                                                    if k6 < ind_homo:
                                                        new_nelec6=new_nelec5+2
                                                    else:
                                                        new_nelec6=new_nelec5
        
                                                    for i7,k7 in enumerate(vir_index):
                                                        if k7 > k6:
                                                            if k7 < ind_homo:
                                                                new_nelec7=new_nelec6+2
                                                            else:
                                                                new_nelec7=new_nelec6
        
                                                            for i8,k8 in enumerate(vir_index):
                                                                if k8 > k7:
                                                                    if k8 < ind_homo:
                                                                        new_nelec8=new_nelec7+2
                                                                    else:
                                                                        new_nelec8=new_nelec7
        
                                                                    new_cas=copy.copy(caslst)
                                                                    new_cas.append(kk)
                                                                    new_cas.append(ll)
                                                                    new_cas.append(lk)
                                                                    new_cas.append(kn)
                                                                    new_cas.append(ko)
                                                                    new_cas.append(k6)
                                                                    new_cas.append(k7)
                                                                    new_cas.append(k8)
                                                                    ncas=len(new_cas)
                                                                    mycas8 = mcscf.CASCI(mf,ncas,new_nelec8)
                                                                    mycas8.max_cycle = 500
                                                                    mycas8.frozen=frozen_lst
        
                                                                    orb_indice1=[i+1 for i in new_cas]
                                                                    mo1 = mycas8.sort_mo(orb_indice1,mo_coeff=mo_coeff)
                                                                    e_tot1=mycas8.kernel(mo1)[0]
                                                                    eight_e_tot[ii][jj][ij][im][io][i6][i7][i8]=e_tot1
                                                                    eight_e_corr[ii][jj][ij][im][io][i6][i7][i8]=e_tot1-e_tot-one_e_corr[ii]-one_e_corr[jj]-one_e_corr[ij]-one_e_corr[im]-one_e_corr[io]-one_e_corr[i6]-one_e_corr[i7]-one_e_corr[i8]-two_e_corr[ii][jj]-two_e_corr[ii][ij]-two_e_corr[ii][im]-two_e_corr[ii][io]-two_e_corr[ii][i6]-two_e_corr[ii][i7]-two_e_corr[ii][i8]-two_e_corr[jj][ij]-two_e_corr[jj][im]-two_e_corr[jj][io]-two_e_corr[jj][i6]-two_e_corr[jj][i7]-two_e_corr[jj][i8]-two_e_corr[ij][im]-two_e_corr[ij][io]-two_e_corr[ij][i6]-two_e_corr[ij][i7]-two_e_corr[ij][i8]-two_e_corr[im][io]-two_e_corr[im][i6]-two_e_corr[im][i7]-two_e_corr[im][i8]-two_e_corr[io][i6]-two_e_corr[io][i7]-two_e_corr[io][i8]-two_e_corr[i6][i7]-two_e_corr[i6][i8]-two_e_corr[i7][i8]-three_e_corr[ii][jj][ij]-three_e_corr[ii][jj][im]-three_e_corr[ii][jj][io]-three_e_corr[ii][jj][i6]-three_e_corr[ii][jj][i7]-three_e_corr[ii][jj][i8]-three_e_corr[ii][ij][im]-three_e_corr[ii][ij][io]-three_e_corr[ii][ij][i6]-three_e_corr[ii][ij][i7]-three_e_corr[ii][ij][i8]-three_e_corr[ii][im][io]-three_e_corr[ii][im][i6]-three_e_corr[ii][im][i7]-three_e_corr[ii][im][i8]-three_e_corr[ii][io][i6]-three_e_corr[ii][io][i7]-three_e_corr[ii][io][i8]-three_e_corr[ii][i6][i7]-three_e_corr[ii][i6][i8]-three_e_corr[ii][i7][i8]-three_e_corr[jj][ij][im]-three_e_corr[jj][ij][io]-three_e_corr[jj][ij][i6]-three_e_corr[jj][ij][i7]-three_e_corr[jj][ij][i8]-three_e_corr[jj][im][io]-three_e_corr[jj][im][i6]-three_e_corr[jj][im][i7]-three_e_corr[jj][im][i8]-three_e_corr[jj][io][i6]-three_e_corr[jj][io][i7]-three_e_corr[jj][io][i8]-three_e_corr[jj][i6][i7]-three_e_corr[jj][i6][i8]-three_e_corr[jj][i7][i8]-three_e_corr[ij][im][io]-three_e_corr[ij][im][i6]-three_e_corr[ij][im][i7]-three_e_corr[ij][im][i8]-three_e_corr[ij][io][i6]-three_e_corr[ij][io][i7]-three_e_corr[ij][io][i8]-three_e_corr[ij][i6][i7]-three_e_corr[ij][i6][i8]-three_e_corr[ij][i7][i8]-three_e_corr[im][io][i6]-three_e_corr[im][io][i7]-three_e_corr[im][io][i8]-three_e_corr[im][i6][i7]-three_e_corr[im][i6][i8]-three_e_corr[im][i7][i8]-three_e_corr[io][i6][i7]-three_e_corr[io][i6][i8]-three_e_corr[io][i7][i8]-three_e_corr[i6][i7][i8]-four_e_corr[ii][jj][ij][im]-four_e_corr[ii][jj][ij][io]-four_e_corr[ii][jj][ij][i6]-four_e_corr[ii][jj][ij][i7]-four_e_corr[ii][jj][ij][i8]-four_e_corr[ii][jj][im][io]-four_e_corr[ii][jj][im][i6]-four_e_corr[ii][jj][im][i7]-four_e_corr[ii][jj][im][i8]-four_e_corr[ii][jj][io][i6]-four_e_corr[ii][jj][io][i7]-four_e_corr[ii][jj][io][i8]-four_e_corr[ii][jj][i6][i7]-four_e_corr[ii][jj][i6][i8]-four_e_corr[ii][jj][i7][i8]-four_e_corr[ii][ij][im][io]-four_e_corr[ii][ij][im][i6]-four_e_corr[ii][ij][im][i7]-four_e_corr[ii][ij][im][i8]-four_e_corr[ii][ij][io][i6]-four_e_corr[ii][ij][io][i7]-four_e_corr[ii][ij][io][i8]-four_e_corr[ii][ij][i6][i7]-four_e_corr[ii][ij][i6][i8]-four_e_corr[ii][ij][i7][i8]-four_e_corr[ii][im][io][i6]-four_e_corr[ii][im][io][i7]-four_e_corr[ii][im][io][i8]-four_e_corr[ii][im][i6][i7]-four_e_corr[ii][im][i6][i8]-four_e_corr[ii][im][i7][i8]-four_e_corr[ii][io][i6][i7]-four_e_corr[ii][io][i6][i8]-four_e_corr[ii][io][i7][i8]-four_e_corr[ii][i6][i7][i8]-four_e_corr[jj][ij][im][io]-four_e_corr[jj][ij][im][i6]-four_e_corr[jj][ij][im][i7]-four_e_corr[jj][ij][im][i8]-four_e_corr[jj][ij][io][i6]-four_e_corr[jj][ij][io][i7]-four_e_corr[jj][ij][io][i8]-four_e_corr[jj][ij][i6][i7]-four_e_corr[jj][ij][i6][i8]-four_e_corr[jj][ij][i7][i8]-four_e_corr[jj][im][io][i6]-four_e_corr[jj][im][io][i7]-four_e_corr[jj][im][io][i8]-four_e_corr[jj][im][i6][i7]-four_e_corr[jj][im][i6][i8]-four_e_corr[jj][im][i7][i8]-four_e_corr[jj][io][i6][i7]-four_e_corr[jj][io][i6][i8]-four_e_corr[jj][io][i7][i8]-four_e_corr[jj][i6][i7][i8]-four_e_corr[ij][im][io][i6]-four_e_corr[ij][im][io][i7]-four_e_corr[ij][im][io][i8]-four_e_corr[ij][im][i6][i7]-four_e_corr[ij][im][i6][i8]-four_e_corr[ij][im][i7][i8]-four_e_corr[ij][io][i6][i7]-four_e_corr[ij][io][i6][i8]-four_e_corr[ij][io][i7][i8]-four_e_corr[ij][i6][i7][i8]-four_e_corr[im][io][i6][i7]-four_e_corr[im][io][i6][i8]-four_e_corr[im][io][i7][i8]-four_e_corr[im][i6][i7][i8]-four_e_corr[io][i6][i7][i8]-five_e_corr[ii][jj][ij][im][io]-five_e_corr[ii][jj][ij][im][i6]-five_e_corr[ii][jj][ij][im][i7]-five_e_corr[ii][jj][ij][im][i8]-five_e_corr[ii][jj][ij][io][i6]-five_e_corr[ii][jj][ij][io][i7]-five_e_corr[ii][jj][ij][io][i8]-five_e_corr[ii][jj][ij][i6][i7]-five_e_corr[ii][jj][ij][i6][i8]-five_e_corr[ii][jj][ij][i7][i8]-five_e_corr[ii][jj][im][io][i6]-five_e_corr[ii][jj][im][io][i7]-five_e_corr[ii][jj][im][io][i8]-five_e_corr[ii][jj][im][i6][i7]-five_e_corr[ii][jj][im][i6][i8]-five_e_corr[ii][jj][im][i7][i8]-five_e_corr[ii][jj][io][i6][i7]-five_e_corr[ii][jj][io][i6][i8]-five_e_corr[ii][jj][io][i7][i8]-five_e_corr[ii][jj][i6][i7][i8]-five_e_corr[ii][ij][im][io][i6]-five_e_corr[ii][ij][im][io][i7]-five_e_corr[ii][ij][im][io][i8]-five_e_corr[ii][ij][im][i6][i7]-five_e_corr[ii][ij][im][i6][i8]-five_e_corr[ii][ij][im][i7][i8]-five_e_corr[ii][ij][io][i6][i7]-five_e_corr[ii][ij][io][i6][i8]-five_e_corr[ii][ij][io][i7][i8]-five_e_corr[ii][ij][i6][i7][i8]-five_e_corr[ii][im][io][i6][i7]-five_e_corr[ii][im][io][i6][i8]-five_e_corr[ii][im][io][i7][i8]-five_e_corr[ii][im][i6][i7][i8]-five_e_corr[ii][io][i6][i7][i8]-five_e_corr[jj][ij][im][io][i6]-five_e_corr[jj][ij][im][io][i7]-five_e_corr[jj][ij][im][io][i8]-five_e_corr[jj][ij][im][i6][i7]-five_e_corr[jj][ij][im][i6][i8]-five_e_corr[jj][ij][im][i7][i8]-five_e_corr[jj][ij][io][i6][i7]-five_e_corr[jj][ij][io][i6][i8]-five_e_corr[jj][ij][io][i7][i8]-five_e_corr[jj][ij][i6][i7][i8]-five_e_corr[jj][im][io][i6][i7]-five_e_corr[jj][im][io][i6][i8]-five_e_corr[jj][im][io][i7][i8]-five_e_corr[jj][im][i6][i7][i8]-five_e_corr[jj][io][i6][i7][i8]-five_e_corr[ij][im][io][i6][i7]-five_e_corr[ij][im][io][i6][i8]-five_e_corr[ij][im][io][i7][i8]-five_e_corr[ij][im][i6][i7][i8]-five_e_corr[ij][io][i6][i7][i8]-five_e_corr[im][io][i6][i7][i8]-six_e_corr[ii][jj][ij][im][io][i6]-six_e_corr[ii][jj][ij][im][io][i7]-six_e_corr[ii][jj][ij][im][io][i8]-six_e_corr[ii][jj][ij][im][i6][i7]-six_e_corr[ii][jj][ij][im][i6][i8]-six_e_corr[ii][jj][ij][im][i7][i8]-six_e_corr[ii][jj][ij][io][i6][i7]-six_e_corr[ii][jj][ij][io][i6][i8]-six_e_corr[ii][jj][ij][io][i7][i8]-six_e_corr[ii][jj][ij][i6][i7][i8]-six_e_corr[ii][jj][im][io][i6][i7]-six_e_corr[ii][jj][im][io][i6][i8]-six_e_corr[ii][jj][im][io][i7][i8]-six_e_corr[ii][jj][im][i6][i7][i8]-six_e_corr[ii][jj][io][i6][i7][i8]-six_e_corr[ii][ij][im][io][i6][i7]-six_e_corr[ii][ij][im][io][i6][i8]-six_e_corr[ii][ij][im][io][i7][i8]-six_e_corr[ii][ij][im][i6][i7][i8]-six_e_corr[ii][ij][io][i6][i7][i8]-six_e_corr[ii][im][io][i6][i7][i8]-six_e_corr[jj][ij][im][io][i6][i7]-six_e_corr[jj][ij][im][io][i6][i8]-six_e_corr[jj][ij][im][io][i7][i8]-six_e_corr[jj][ij][im][i6][i7][i8]-six_e_corr[jj][ij][io][i6][i7][i8]-six_e_corr[jj][im][io][i6][i7][i8]-six_e_corr[ij][im][io][i6][i7][i8]-seven_e_corr[ii][jj][ij][im][io][i6][i7]-seven_e_corr[ii][jj][ij][im][io][i6][i8]-seven_e_corr[ii][jj][ij][im][io][i7][i8]-seven_e_corr[ii][jj][ij][im][i6][i7][i8]-seven_e_corr[ii][jj][ij][io][i6][i7][i8]-seven_e_corr[ii][jj][im][io][i6][i7][i8]-seven_e_corr[ii][ij][im][io][i6][i7][i8]-seven_e_corr[jj][ij][im][io][i6][i7][i8]
                                                                    with open('eight-body.txt','a') as f:
                                                                        np.savetxt(f,[[kk,ll,lk,kn,ko,k6,k7,k8,e_tot1,eight_e_corr[ii][jj][ij][im][io][i6][i7][i8]]],fmt='%d,%d,%d,%d,%d,%d,%d,%d,%f,%f')
                                                                    del new_cas
        
        print('eight-all :', base_e_corr+np.sum(one_e_corr)+np.sum(two_e_corr)+np.sum(three_e_corr)+np.sum(four_e_corr)+np.sum(five_e_corr)+np.sum(six_e_corr)+np.sum(seven_e_corr)+np.sum(eight_e_corr))
    sys.stdout.flush()


