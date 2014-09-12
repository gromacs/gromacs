/*
 * $Id: gentop.c,v 1.26 2009/05/20 10:48:03 spoel Exp $
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
/* gaussian_integrals.c (c) 2010 Paul J. van Maaren and David van der Spoel */
#include "gmxpre.h"
#include <stdio.h>
#include <math.h>
#include "coulombintegrals.h"

static double sqr(double x)
{
    return x*x;
}

double Nuclear_GG(double r, double zeta)
{
    /* This routine may be called with zeta 0.
     * In that case it is a simple 1/r interaction.
     */
    if (zeta == 0)
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return 1.0/r;
        }
    }
    else if (r == 0)
    {
        return 2.0*zeta/sqrt(M_PI);
    }
    else
    {
        return erf(zeta*r)/r;
    }
}

double Coulomb_GG(double r, double zi, double zj)
{
    double zeff;

    /* This routine may be called with one or both zeta 0.
     * In that case it is either a Nuclear_GG or simple 1/r interaction.
     */
    if ((zi == 0) && (zj == 0))
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return 1/r;
        }
    }
    else
    {
        if (zi == 0)
        {
            zeff = zj;
        }
        else if (zj == 0)
        {
            zeff = zi;
        }
        else
        {
            zeff = zi*zj/sqrt(sqr(zi)+sqr(zj));
        }

        return Nuclear_GG(r, zeff);
    }
}

double DNuclear_GG(double r, double z)
{
    double rz = r * z;
    double dngg;

    if (z == 0)
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1.0/sqr(r);
        }
    }
    else if (r == 0)
    {
        dngg = 0;
    }
    else
    {
        dngg = (2.0/sqrt(M_PI))*exp(-sqr(rz))*(z/r) - erf(rz)/sqr(r);
    }

    return dngg;
}

double DCoulomb_GG(double r, double zi, double zj)
{
    double zeff;

    if ((zi == 0) && (zj == 0))
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1/sqr(r);
        }
    }
    else
    {
        if (zi == 0)
        {
            zeff = zj;
        }
        else if (zj == 0)
        {
            zeff = zi;
        }
        else
        {
            zeff = zi*zj/sqrt(sqr(zi)+sqr(zj));
        }

        return DNuclear_GG(r, zeff);
    }
}
