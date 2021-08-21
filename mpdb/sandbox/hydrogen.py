"""
Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import CoolProp as cP

c1="Hydrogen"
c2="Nitrogen"
CAS_c1=cP.CoolProp.get_fluid_param_string(c1,'CAS')
CAS_c2=cP.CoolProp.get_fluid_param_string(c2,'CAS')
print("betaT: %.4f;" % float(cP.CoolProp.get_mixture_binary_pair_data(CAS_c1,CAS_c2,'betaT')))
print("betaV: %.4f;" % float(cP.CoolProp.get_mixture_binary_pair_data(CAS_c1,CAS_c2,'betaV')))

#https://github.com/ibell/coolprop/blob/master/dev/mixtures/KunzWagner2012_TableA8.txt
#cP.CoolProp.set_mixture_binary_pair_data(CAS_c1,CAS_c2,'betaT',0.972532065)
#cP.CoolProp.set_mixture_binary_pair_data(CAS_c1,CAS_c2,'betaV',0.946134337)

"""
bip = [["Nitrogen",,],
       ["Oxygen",,],
       ["Argon",,],
       ["Water",,],
       ["CarbonDioxide",,],
       ["CarbonMonoxide",,],
       ["Hydrogen","Methane",0.0705],
       ["Hydrogen","Nitrogen",0.103],
       ["Methane",,],
       ["Ammonia",,],
       ["HydrogenSulfide",,]]
 """
"""
c1="Hydrogen"
c2="Methane"
dataHeader=["T /K", "P /bar", "%s lMF" % c1 , "%s vMF" % c1]
data=[[103.15,19.961,0.0142, 0.9530],
#      [103.15,41.34,0.0337,0.9670],
#      [103.15,61.2,0.0569,0.9690],
#      [103.15,81.71,0.0708,0.9700],
#      [103.15,91.192,0.0819,0.9690],
      [123.15,10.234,0.0073,0.6830],
      [123.15,20.265,0.0192,0.8180],
      [123.15,40.327,0.0385,0.8860],
      [123.15,60.592,0.0616,0.9040],
      [123.15,81.769,0.0816,0.9110],
      [123.15,101.831,0.1010,0.9080],
      [143.05,10.740,0.0031,0.2140],
      [143.05,20.164,0.0168,0.5410],
      [143.05,40.732,0.0477,0.7210],
      [143.05,59.883,0.0684,0.7790],
      [143.05,79.742,0.0992,0.8020],
      [143.05,101.527,0.1270,0.8100],
      [172.05,87.544,0.1360,0.4350],
      [172.05,103.553,0.2010,0.4290],
      [173.65,35.666,0.0168,0.1480],
      [173.65,40.631,0.0257,0.2060],
      [173.65,60.693,0.0709,0.3620],
      [173.65,71.028,0.0938,0.4060],
      [173.65,81.667,0.124,0.4240],
      [173.65,83.795,0.1310,0.4240],
      [173.65,100.514,0.1900,0.4200],
      [173.65,103.351,0.2050,0.4110],
      [173.65,108.316,0.2250,0.4010]]
#print(tabulate(data,headers=dataHeader)+'\n')
"""
"""
c1="Hydrogen"
c2="Nitrogen"
dataHeader=["T /K", "P /bar", "%s lMF" % c1 , "%s vMF" % c1]
data = [
    [90.00,3.595,0.0, 0.0],
    [89.99,8.436,0.0119,0.5217],
    [90.01,10.066,0.0159,0.5855],
    [90.04,10.209,0.0164,0.5903],
    [90.03,11.853,0.0207,0.6367],
    [90.01,12.976,0.0234,0.6619],
    [90.02,13.566,0.0248,0.6732],
    [90.06,24.420,0.0534,0.7800],
    [90.04,45.954,0.1115,0.8304],
    [95.00,5.408,0.0,0.0],
    [95.02,11.091,0.0149,0.4492],
    [95.00,26.216,0.0566,0.7005],
    [94.99,45.328,0.1121,0.7643]
]
#print(tabulate(data,headers=dataHeader)+'\n')
"""
