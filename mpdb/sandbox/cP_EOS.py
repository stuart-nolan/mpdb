"""
cP_EOS.py

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import pdb
import CoolProp.CoolProp as cP
from CoolProp import AbstractState as cPAS
from tabulate import tabulate
from scipy.optimize import minimize

def cEOS_fit_kij(kij, data, cPAS_EOS):
    cPAS_EOS.set_binary_interaction_double(0,1,"kij",kij[0])
    (DPpP,Dy,outData) = deltaVar(data,cPAS_EOS)
    return DPpP # fit based on DPpP

def deltaVar(inData, cPAS_EOS):
    DPpP = 0 # sum_over_i( abs(P_exp_i- P_EOS_i)/P_exp_i )
    Dy = 0 # sum_over_i( abs(vMF_exp_c1_i - vMF_EOS_c1_i )
    outData = []
    for row in inData:
        T_exp = row[0]
        P_exp = row[1]
        lMF_exp_c1 = row[2] # c1 liquid Mole Fraction
        vMF_exp_c1 = row[3] # c1 vapor Mole Fraction

        cPAS_EOS.set_mole_fractions([lMF_exp_c1, 1-lMF_exp_c1]) 
        #if cPAS_cEOS.phase() in [0,6]: # 0 is liquid phase
        cPAS_EOS.update(cP.QT_INPUTS, 0, T_exp) # VLE, 0 = 100% liquid phase
        P_EOS = cPAS_EOS.p()/100000 # bar
        DPpP = DPpP + abs(P_exp - P_EOS)/P_exp
        vMFs_EOS = cPAS_EOS.mole_fractions_vapor() # EOS vMFs = [c1, c2]
        Dy = Dy + abs(vMF_exp_c1 - vMFs_EOS[0])
        outData.append([T_exp, lMF_exp_c1, P_exp,
                        "{0:.3f}".format(P_EOS),
                        "{0:.3f}".format(vMF_exp_c1),
                        "{0:.3f}".format(vMFs_EOS[0])])
        
    DPpP = DPpP*100/len(data)
    Dy = Dy*100/len(data)
    return (DPpP, Dy, outData)

"""
from CoolProp/include/DataStructures.h

0 iphase_liquid, < Subcritical liquid
1 iphase_supercritical, < Supercritical (p > pc, T > Tc)
2 iphase_supercritical_gas, < Supercritical gas (p < pc, T > Tc)
3 iphase_supercritical_liquid, < Supercritical liquid (p > pc, T < Tc)
4 iphase_critical_point, < At the critical point
5 iphase_gas, < Subcritical gas
6 iphase_twophase, < Twophase
7 iphase_unknown, < Unknown phase
8 iphase_not_imposed

import CoolProp
CoolProp.iphase_twophase
Out: 6
"""

cPFluids = dict([(fluid, cP.get_fluid_param_string(fluid,"CAS")) for fluid in 
                 cP.get_global_param_string("fluids_list").split(',')])

# H2, N2, CO ternary data, Table VI, pA-4 (p78), Eubanks, 1957
# REF: https://scholarship.rice.edu/bitstream/handle/1911/18254/3079688.PDF?sequence=1&isAllowed=y
c1="Hydrogen"
c2="Nitrogen"
rawdata_header=["T /degF", "P /psia", "%s lMF" % c1 , "%s vMF" % c1]
rawdata = [
    [-310,315,0.0487, 0.8655],
    [-310,500,0.0763, 0.8948],
#    [-310,1400,0.2488, 0.8622],
#    [-310,2000,0.3446, 0.7977],
    [-280,315,0.0377, 0.5509],
    [-280,500,0.0741, 0.6686],
    [-280,800,0.1384, 0.7205],
    [-280,1100,0.2092, 0.7070],
    [-280,1400,0.3221, 0.6462],
]
data_header=["T /K", "P /bar", "%s lMF" % c1 , "%s vMF" % c1]
data = []
psia2bar = lambda P : P/14.503774
degF2K = lambda T : 5/9*(T - 32) + 273.15

for row in rawdata:
    data.append([degF2K(row[0]),psia2bar(row[1]),row[2],row[3]])

#print(tabulate(data,headers=data_header)+'\n')

EOS="PR"
cPAS_cEOS = cPAS(EOS, c1+"&"+c2)
res = minimize(cEOS_fit_kij, 0.1, bounds=[(-0.2,0.5)], args=(data, cPAS_cEOS))
kij=res.x[0]

#kij=0.0864 # EOS="PR"
#kij=0.0641 # EOS="SRK"

cPAS_cEOS.set_binary_interaction_double(0,1,"kij",kij)
(DPpP,Dy,outData) = deltaVar(data,cPAS_cEOS)

outData_header=["T_exp /K", "c1 lMF_exp", "P_exp /bar", "P_%s /bar" % EOS,
                "c1 vMF_exp", "c1 vMF_%s" % EOS]
print("\nc1: %s;" % c1)
print("c2: %s;" % c2)
print("kij: %.4f; DPpP_%s: %.2f; Dy_%s: %.2f;\n" % (kij, EOS, DPpP, EOS, Dy))
print(tabulate(outData,headers=outData_header)+'\n')

CAS_c1=cPFluids[c1]
CAS_c2=cPFluids[c2]
# [c2-c1 betaT gammaT betaV gammaV], Table A8, ref 31, Kunz & Wagner 2012
# REF: https://github.com/CoolProp/CoolProp/blob/master/dev/mixtures/KunzWagner2012_TableA8.txt
#cP.set_mixture_binary_pair_data(CAS_c1,CAS_c2,'betaT',0.972532065)
#cP.set_mixture_binary_pair_data(CAS_c1,CAS_c2,'betaV',0.946134337)

EOS = "HEOS"
cPAS_HEOS = cP.AbstractState(EOS,c1+"&"+c2)
"""
(DPpP,Dy,outData) = deltaVar(data,cPAS_HEOS)
# pdb.pm() from ipython...
outData_header=["T_exp /K", "c1 lMF_exp", "P_exp /bar", "P_%s /bar" % EOS,
                "c1 vMF_exp", "c1 vMF_%s" % EOS]
print("\nc1: %s;" % c1)
print("c2: %s;" % c2)
print("DPpP_%s: %.2f; Dy_%s: %.2f;\n" % (EOS, DPpP, EOS, Dy))
print(tabulate(outData,headers=outData_header)+'\n')
"""
