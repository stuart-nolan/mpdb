"""
Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.

    heatOfRxn(dictionary comp, dictionary db, float T):
            
    Parameters:
        comp, dictionary of components keys
              values: dictionary {stoic: [+,-] float,
                                  state: [v,l,s]}
                      stoic key: negative values for reactants
                      state: v = vapor, l = liquid, s = solid

        db, database of pure component material properites, each comp key must
            be a key in db
        T, temperature in K
    
    Returns:
        heat of reaction at T in J/mol
"""

# check for free energy formation for each comp in db
# check for vapor, liquid, solid heat capacity data for each comp in db
# check for state changes between reference temp 298.15 K and T
