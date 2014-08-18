/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "units.h"

#include <stdio.h>

#include "gromacs/utility/cstringutil.h"

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
