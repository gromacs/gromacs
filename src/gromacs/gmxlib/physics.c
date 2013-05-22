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

double convert2gmx(double x, int unit)
{
    switch (unit)
    {
        case eg2cAngstrom:
            return x*A2NM;
        case eg2cNm:
            return x;
        case eg2cBohr:
            return x*BOHR2NM;
        case eg2cKcal_Mole:
            return x/CAL2JOULE;
        case eg2cHartree:
            return x*ONE_4PI_EPS0/BOHR2NM;
        case eg2cHartree_e:
            return x*ONE_4PI_EPS0/BOHR2NM;
        case eg2cAngstrom3:
            return x*A2NM*A2NM*A2NM;
        case eg2cCoulomb:
            return x/E_CHARGE;
        case eg2cDebye:
            return x*DEBYE2ENM;
        case eg2cElectron:
            return x;
        case eg2cBuckingham:
            return x*A2NM*DEBYE2ENM;
        default:
            fprintf(stderr, "Unknown unit %d, not converting.\n", unit);
    }
    return x;
}

double gmx2convert(double x, int unit)
{
    switch (unit)
    {
        case eg2cAngstrom:
            return x/A2NM;
        case eg2cNm:
            return x;
        case eg2cBohr:
            return x/BOHR2NM;
        case eg2cKcal_Mole:
            return x*CAL2JOULE;
        case eg2cHartree:
            return x/(ONE_4PI_EPS0/BOHR2NM);
        case eg2cHartree_e:
            return x/(ONE_4PI_EPS0/BOHR2NM);
        case eg2cAngstrom3:
            return x/(A2NM*A2NM*A2NM);
        case eg2cCoulomb:
            return x*E_CHARGE;
        case eg2cDebye:
            return x/DEBYE2ENM;
        case eg2cElectron:
            return x;
        case eg2cBuckingham:
            return x/(A2NM*DEBYE2ENM);
        default:
            fprintf(stderr, "Unknown unit %d, not converting.\n", unit);
    }
    return x;
}

/* This has to have the same order as the enums. */
static const char *eg2c_names[eg2cNR] = {
    "Angstrom", "Nm", "Bohr", "Kcal_Mole",
    "Hartree", "Hartree_e", "Angstrom3", "Coulomb",
    "Debye", "Electron", "Buckingham"
};

int string2unit(char *string)
{
    int i;

    for (i = 0; (i < eg2cNR); i++)
    {
        if (gmx_strcasecmp(string, eg2c_names[i]) == 0)
        {
            return i;
        }
    }
    return -1;
}

const char *unit2string(int unit)
{
    if ((unit >= 0) && (unit < eg2cNR))
    {
        return eg2c_names[unit];
    }

    return NULL;
}
