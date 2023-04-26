#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:10:01 2023

@authors:    
Ryan Boukrouche (RB)
"""

import numpy as np
import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles

#Moist adjustment routine.
def moist_adj(Tmid,pmid,pint,conv_timescale):
    dT_conv = np.zeros(len(pmid))
    Tmid_cc = np.zeros(len(pmid))
    for i in range(len(pmid)):
        Tmid_cc[i] = Tmid[i]
    nb_convsteps = 10
    for i_conv in range(nb_convsteps):
        did_adj = False
        #------------------- L = cst -------------------
        for i in range(len(Tmid_cc)-1): #Downward pass
            if (Tmid_cc[i] < ga.Tdew('H2O',pmid[i])):
                Tmid_cc[i]=ga.Tdew('H2O',pmid[i])
                did_adj = True
        for i in range(len(Tmid_cc)-2,-1,-1): #Upward pass
            if (Tmid_cc[i] < ga.Tdew('H2O',pmid[i])):
                Tmid_cc[i]=ga.Tdew('H2O',pmid[i])
                did_adj = True

        # If no adjustment required, exit the loop
        if (did_adj == False):
            break

    # Change in temperature is Tmid_cc - Tmid
    dT_conv[:] = (Tmid_cc[:] - Tmid[:])/conv_timescale
    return dT_conv
