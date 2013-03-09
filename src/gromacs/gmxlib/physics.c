/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop_util.c,v 1.51 2009/06/01 06:13:18 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#include <stdio.h>
#include "string2.h"
#include "physics.h"
	
double convert2gmx(double x,int unit)
{
    switch (unit) 
    {
    case eg2c_Angstrom:
        return x*A2NM;
    case eg2c_nm:
        return x;
    case eg2c_pm:
        return 0.001*x;
    case eg2c_Bohr:
        return x*BOHR2NM;
    case eg2c_kcal_mole:
        return x*CAL2JOULE;
    case eg2c_kJ_mole:
        return x;
    case eg2c_Hartree:
        return x*ONE_4PI_EPS0/BOHR2NM;
    case eg2c_Hartree_e:
        return x*ONE_4PI_EPS0/BOHR2NM;
    case eg2c_Angstrom3:
        return x*A2NM*A2NM*A2NM;
    case eg2c_Coulomb:
        return x/E_CHARGE;
    case eg2c_Debye:
        return x*DEBYE2ENM;
    case eg2c_Electron:
        return x;
    case eg2c_Buckingham:
        return x*A2NM*DEBYE2ENM;
    default:
        fprintf(stderr,"Unknown unit %d, not converting.\n",unit);
    }  
    return x;
}

double gmx2convert(double x,int unit)
{
    switch (unit) 
    {
    case eg2c_Angstrom:
        return x/A2NM;
    case eg2c_nm:
        return x;
    case eg2c_pm:
        return 1000*x;
    case eg2c_Bohr:
        return x/BOHR2NM;
    case eg2c_kcal_mole:
        return x/CAL2JOULE;
    case eg2c_kJ_mole:
        return x;
    case eg2c_Hartree:
        return x/(ONE_4PI_EPS0/BOHR2NM);
    case eg2c_Hartree_e:
        return x/(ONE_4PI_EPS0/BOHR2NM);
    case eg2c_Angstrom3:
        return x/(A2NM*A2NM*A2NM);
    case eg2c_Coulomb:
        return x*E_CHARGE;
    case eg2c_Debye:
        return x/DEBYE2ENM;
    case eg2c_Electron:
        return x;
    case eg2c_Buckingham:
        return x/(A2NM*DEBYE2ENM);
    default:
        fprintf(stderr,"Unknown unit %d, not converting.\n",unit);
    }  
    return x;
}

/* This has to have the same order as the enums. */
static const char *eg2c_names[eg2c_NR] = {
    "Angstrom", "nm", "pm", "Bohr", "kcal/mol", "kJ/mol",
    "Hartree", "Hartree/e", "Angstrom3", "Coulomb",
    "Debye", "Electron", "Buckingham" 
};

int string2unit(const char *string)
{
    int i;
    
    for(i=0; (i<eg2c_NR); i++)
        if (gmx_strcasecmp(string,eg2c_names[i]) == 0)
            return i;
    return -1;
}

const char *unit2string(int unit)
{
    if ((unit >= 0) && (unit < eg2c_NR))
        return eg2c_names[unit];
        
    return NULL;
}
