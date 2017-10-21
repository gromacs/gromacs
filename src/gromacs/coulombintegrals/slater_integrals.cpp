/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
/* slater_integrals.cpp (c) 2010 Paul J. van Maaren and David van der Spoel */
#include "gmxpre.h"

#include <stdio.h>
#include <iostream>
#include <cmath>

#include "coulombintegrals.h"
#include "slater_low.h"

#ifdef HAVE_LIBCLN

cl_R Nuclear_1S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = 1LL/r - (1LL + r*xi)/(exp(2LL*r*xi)*r)

    ;

    return S;
}

cl_R Nuclear_2S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = 1LL/r - (6LL + 9LL*r*xi + 6LL*Power(r, 2LL)*Power(xi, 2LL) + 2LL*Power(r, 3LL)*Power(xi, 3LL))/

        (6LL*exp(2LL*r*xi)*r)

    ;

    return S;
}

cl_R Nuclear_3S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = 1LL/r - (45LL + 75LL*r*xi + 60LL*Power(r, 2LL)*Power(xi, 2LL) +

                 30LL*Power(r, 3LL)*Power(xi, 3LL) + 10LL*Power(r, 4LL)*Power(xi, 4LL) +

                 2LL*Power(r, 5LL)*Power(xi, 5LL))/(45LL*exp(2LL*r*xi)*r)

    ;

    return S;
}

cl_R Nuclear_4S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = 1LL/r - (1260LL + 2205LL*r*xi + 1890LL*Power(r, 2LL)*Power(xi, 2LL) +

                 1050LL*Power(r, 3LL)*Power(xi, 3LL) + 420LL*Power(r, 4LL)*Power(xi, 4LL) +

                 126LL*Power(r, 5LL)*Power(xi, 5LL) + 28LL*Power(r, 6LL)*Power(xi, 6LL) +

                 4LL*Power(r, 7LL)*Power(xi, 7LL))/(1260LL*exp(2LL*r*xi)*r)

    ;

    return S;
}

cl_R Nuclear_5S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = 1LL/r - (14175LL + 25515LL*r*xi + 22680LL*Power(r, 2LL)*Power(xi, 2LL) +

                 13230LL*Power(r, 3LL)*Power(xi, 3LL) + 5670LL*Power(r, 4LL)*Power(xi, 4LL) +

                 1890LL*Power(r, 5LL)*Power(xi, 5LL) + 504LL*Power(r, 6LL)*Power(xi, 6LL) +

                 108LL*Power(r, 7LL)*Power(xi, 7LL) + 18LL*Power(r, 8LL)*Power(xi, 8LL) +

                 2LL*Power(r, 9LL)*Power(xi, 9LL))/(14175LL*exp(2LL*r*xi)*r)

    ;

    return S;
}

cl_R Nuclear_6S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = 1LL/r - (935550LL + 1715175LL*r*xi + 1559250LL*Power(r, 2LL)*Power(xi, 2LL) +

                 935550LL*Power(r, 3LL)*Power(xi, 3LL) + 415800LL*Power(r, 4LL)*Power(xi, 4LL) +

                 145530LL*Power(r, 5LL)*Power(xi, 5LL) + 41580LL*Power(r, 6LL)*Power(xi, 6LL) +

                 9900LL*Power(r, 7LL)*Power(xi, 7LL) + 1980LL*Power(r, 8LL)*Power(xi, 8LL) +

                 330LL*Power(r, 9LL)*Power(xi, 9LL) + 44LL*Power(r, 10LL)*Power(xi, 10LL) +

                 4LL*Power(r, 11LL)*Power(xi, 11LL))/(935550LL*exp(2LL*r*xi)*r)

    ;

    return S;
}

cl_R DNuclear_1S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = Power(r, -2LL) - (1LL + 2LL*r*xi + 2LL*Power(r, 2LL)*Power(xi, 2LL))/

        (exp(2LL*r*xi)*Power(r, 2LL))

    ;

    return S;
}

cl_R DNuclear_2S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = Power(r, -2LL) - (3LL + 6LL*r*xi + 6LL*Power(r, 2LL)*Power(xi, 2LL) +

                          4LL*Power(r, 3LL)*Power(xi, 3LL) + 2LL*Power(r, 4LL)*Power(xi, 4LL))/

        (3LL*exp(2LL*r*xi)*Power(r, 2LL))

    ;

    return S;
}

cl_R DNuclear_3S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = Power(r, -2LL) - (45LL + 90LL*r*xi + 90LL*Power(r, 2LL)*Power(xi, 2LL) +

                          60LL*Power(r, 3LL)*Power(xi, 3LL) + 30LL*Power(r, 4LL)*Power(xi, 4LL) +

                          12LL*Power(r, 5LL)*Power(xi, 5LL) + 4LL*Power(r, 6LL)*Power(xi, 6LL))/

        (45LL*exp(2LL*r*xi)*Power(r, 2LL))

    ;

    return S;
}

cl_R DNuclear_4S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = Power(r, -2LL) - (315LL + 630LL*r*xi + 630LL*Power(r, 2LL)*Power(xi, 2LL) +

                          420LL*Power(r, 3LL)*Power(xi, 3LL) + 210LL*Power(r, 4LL)*Power(xi, 4LL) +

                          84LL*Power(r, 5LL)*Power(xi, 5LL) + 28LL*Power(r, 6LL)*Power(xi, 6LL) +

                          8LL*Power(r, 7LL)*Power(xi, 7LL) + 2LL*Power(r, 8LL)*Power(xi, 8LL))/

        (315LL*exp(2LL*r*xi)*Power(r, 2LL))

    ;

    return S;
}

cl_R DNuclear_5S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = Power(r, -2LL) - (14175LL + 28350LL*r*xi + 28350LL*Power(r, 2LL)*Power(xi, 2LL) +

                          18900LL*Power(r, 3LL)*Power(xi, 3LL) + 9450LL*Power(r, 4LL)*Power(xi, 4LL) +

                          3780LL*Power(r, 5LL)*Power(xi, 5LL) + 1260LL*Power(r, 6LL)*Power(xi, 6LL) +

                          360LL*Power(r, 7LL)*Power(xi, 7LL) + 90LL*Power(r, 8LL)*Power(xi, 8LL) +

                          20LL*Power(r, 9LL)*Power(xi, 9LL) + 4LL*Power(r, 10LL)*Power(xi, 10LL))/

        (14175LL*exp(2LL*r*xi)*Power(r, 2LL))

    ;

    return S;
}

cl_R DNuclear_6S(cl_R r, cl_R xi)
{
    cl_R S = ZERO;
    S = Power(r, -2LL) - (467775LL + 935550LL*r*xi + 935550LL*Power(r, 2LL)*Power(xi, 2LL) +

                          623700LL*Power(r, 3LL)*Power(xi, 3LL) + 311850LL*Power(r, 4LL)*Power(xi, 4LL) +

                          124740LL*Power(r, 5LL)*Power(xi, 5LL) + 41580LL*Power(r, 6LL)*Power(xi, 6LL) +

                          11880LL*Power(r, 7LL)*Power(xi, 7LL) + 2970LL*Power(r, 8LL)*Power(xi, 8LL) +

                          660LL*Power(r, 9LL)*Power(xi, 9LL) + 132LL*Power(r, 10LL)*Power(xi, 10LL) +

                          24LL*Power(r, 11LL)*Power(xi, 11LL) + 4LL*Power(r, 12LL)*Power(xi, 12LL))/

        (467775LL*exp(2LL*r*xi)*Power(r, 2LL))

    ;

    return S;
}

typedef cl_R t_slater_SS_func (cl_R r, cl_R xi, cl_R xj);
typedef cl_R t_slater_NS_func (cl_R r, cl_R xi);
t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]) = {
    {  Slater_1S_1S,  Slater_1S_2S,  Slater_1S_3S ,  Slater_1S_4S,  Slater_1S_5S,  Slater_1S_6S},
    {  Slater_2S_1S,  Slater_2S_2S,  Slater_2S_3S ,  Slater_2S_4S,  Slater_2S_5S,  Slater_2S_6S},
    {  Slater_3S_1S,  Slater_3S_2S,  Slater_3S_3S ,  Slater_3S_4S,  Slater_3S_5S,  Slater_3S_6S},
    {  Slater_4S_1S,  Slater_4S_2S,  Slater_4S_3S,   Slater_4S_4S,  Slater_4S_5S,  Slater_4S_6S},
    {  Slater_5S_1S,  Slater_5S_2S,  Slater_5S_3S,   Slater_5S_4S,  Slater_5S_5S,  Slater_5S_6S},
    {  Slater_6S_1S,  Slater_6S_2S,  Slater_6S_3S,   Slater_6S_4S,  Slater_6S_5S,  Slater_6S_6S}
};

t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]) = {
    {  DSlater_1S_1S,  DSlater_1S_2S,  DSlater_1S_3S ,  DSlater_1S_4S,  DSlater_1S_5S,  DSlater_1S_6S},
    {  DSlater_2S_1S,  DSlater_2S_2S,  DSlater_2S_3S ,  DSlater_2S_4S,  DSlater_2S_5S,  DSlater_2S_6S},
    {  DSlater_3S_1S,  DSlater_3S_2S,  DSlater_3S_3S ,  DSlater_3S_4S,  DSlater_3S_5S,  DSlater_3S_6S},
    {  DSlater_4S_1S,  DSlater_4S_2S,  DSlater_4S_3S,   DSlater_4S_4S,  DSlater_4S_5S,  DSlater_4S_6S},
    {  DSlater_5S_1S,  DSlater_5S_2S,  DSlater_5S_3S,   DSlater_5S_4S,  DSlater_5S_5S,  DSlater_5S_6S},
    {  DSlater_6S_1S,  DSlater_6S_2S,  DSlater_6S_3S,   DSlater_6S_4S,  DSlater_6S_5S,  DSlater_6S_6S}
};

t_slater_NS_func (*Slater_NS[SLATER_MAX]) = {
    Nuclear_1S,  Nuclear_2S,  Nuclear_3S ,  Nuclear_4S,  Nuclear_5S,  Nuclear_6S
};

t_slater_NS_func (*DSlater_NS[SLATER_MAX]) = {
    DNuclear_1S,  DNuclear_2S,  DNuclear_3S ,  DNuclear_4S,  DNuclear_5S,  DNuclear_6S
};

static char *my_ftoa(double d)
{
    static char buf[256];
    sprintf(buf, "%f", d);
    if (strchr(buf, '.') == nullptr)
    {
        strcat(buf, ".0");
    }
    strcat(buf, "_80");
    return buf;
}

#else
/* NOT HAVE_LIBCLN */

double Nuclear_1S(double r, double xi)
{
    double S = 0;
    S = 1/r - (1 + r*xi)/(exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_2S(double r, double xi)
{
    double S = 0;
    S = 1/r - (6 + 9*r*xi + 6*power(r, 2)*power(xi, 2) + 2*power(r, 3)*power(xi, 3))/

        (6*exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_3S(double r, double xi)
{
    double S = 0;
    S = 1/r - (45 + 75*r*xi + 60*power(r, 2)*power(xi, 2) +

                 30*power(r, 3)*power(xi, 3) + 10*power(r, 4)*power(xi, 4) +

                 2*power(r, 5)*power(xi, 5))/(45*exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_4S(double r, double xi)
{
    double S = 0;
    S = 1/r - (1260 + 2205*r*xi + 1890*power(r, 2)*power(xi, 2) +

                 1050*power(r, 3)*power(xi, 3) + 420*power(r, 4)*power(xi, 4) +

                 126*power(r, 5)*power(xi, 5) + 28*power(r, 6)*power(xi, 6) +

                 4*power(r, 7)*power(xi, 7))/(1260*exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_5S(double r, double xi)
{
    double S = 0;
    S = 1/r - (14175 + 25515*r*xi + 22680*power(r, 2)*power(xi, 2) +

                 13230*power(r, 3)*power(xi, 3) + 5670*power(r, 4)*power(xi, 4) +

                 1890*power(r, 5)*power(xi, 5) + 504*power(r, 6)*power(xi, 6) +

                 108*power(r, 7)*power(xi, 7) + 18*power(r, 8)*power(xi, 8) +

                 2*power(r, 9)*power(xi, 9))/(14175*exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_6S(double r, double xi)
{
    double S = 0;
    S = 1/r - (935550 + 1715175*r*xi + 1559250*power(r, 2)*power(xi, 2) +

                 935550*power(r, 3)*power(xi, 3) + 415800*power(r, 4)*power(xi, 4) +

                 145530*power(r, 5)*power(xi, 5) + 41580*power(r, 6)*power(xi, 6) +

                 9900*power(r, 7)*power(xi, 7) + 1980*power(r, 8)*power(xi, 8) +

                 330*power(r, 9)*power(xi, 9) + 44*power(r, 10)*power(xi, 10) +

                 4*power(r, 11)*power(xi, 11))/(935550*exp(2*r*xi)*r)

    ;

    return S;
}

double DNuclear_1S(double r, double xi)
{
    double S = 0;
    S = power(r, -2) - (1 + 2*r*xi + 2*power(r, 2)*power(xi, 2))/

        (exp(2*r*xi)*power(r, 2))

    ;

    return S;
}

double DNuclear_2S(double r, double xi)
{
    double S = 0;
    S = power(r, -2) - (3 + 6*r*xi + 6*power(r, 2)*power(xi, 2) +

                          4*power(r, 3)*power(xi, 3) + 2*power(r, 4)*power(xi, 4))/

        (3*exp(2*r*xi)*power(r, 2))

    ;

    return S;
}

double DNuclear_3S(double r, double xi)
{
    double S = 0;
    S = power(r, -2) - (45 + 90*r*xi + 90*power(r, 2)*power(xi, 2) +

                          60*power(r, 3)*power(xi, 3) + 30*power(r, 4)*power(xi, 4) +

                          12*power(r, 5)*power(xi, 5) + 4*power(r, 6)*power(xi, 6))/

        (45*exp(2*r*xi)*power(r, 2))

    ;

    return S;
}

double DNuclear_4S(double r, double xi)
{
    double S = 0;
    S = power(r, -2) - (315 + 630*r*xi + 630*power(r, 2)*power(xi, 2) +

                          420*power(r, 3)*power(xi, 3) + 210*power(r, 4)*power(xi, 4) +

                          84*power(r, 5)*power(xi, 5) + 28*power(r, 6)*power(xi, 6) +

                          8*power(r, 7)*power(xi, 7) + 2*power(r, 8)*power(xi, 8))/

        (315*exp(2*r*xi)*power(r, 2))

    ;

    return S;
}

double DNuclear_5S(double r, double xi)
{
    double S = 0;
    S = power(r, -2) - (14175 + 28350*r*xi + 28350*power(r, 2)*power(xi, 2) +

                          18900*power(r, 3)*power(xi, 3) + 9450*power(r, 4)*power(xi, 4) +

                          3780*power(r, 5)*power(xi, 5) + 1260*power(r, 6)*power(xi, 6) +

                          360*power(r, 7)*power(xi, 7) + 90*power(r, 8)*power(xi, 8) +

                          20*power(r, 9)*power(xi, 9) + 4*power(r, 10)*power(xi, 10))/

        (14175*exp(2*r*xi)*power(r, 2))

    ;

    return S;
}

double DNuclear_6S(double r, double xi)
{
    double S = 0;
    S = power(r, -2) - (467775 + 935550*r*xi + 935550*power(r, 2)*power(xi, 2) +

                          623700*power(r, 3)*power(xi, 3) + 311850*power(r, 4)*power(xi, 4) +

                          124740*power(r, 5)*power(xi, 5) + 41580*power(r, 6)*power(xi, 6) +

                          11880*power(r, 7)*power(xi, 7) + 2970*power(r, 8)*power(xi, 8) +

                          660*power(r, 9)*power(xi, 9) + 132*power(r, 10)*power(xi, 10) +

                          24*power(r, 11)*power(xi, 11) + 4*power(r, 12)*power(xi, 12))/

        (467775*exp(2*r*xi)*power(r, 2))

    ;

    return S;
}

typedef double t_slater_SS_func (double r, double xi, double xj);
typedef double t_slater_NS_func (double r, double xi);
t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]) = {
    {  Slater_1S_1S,  Slater_1S_2S,  Slater_1S_3S},
    {  Slater_2S_1S,  Slater_2S_2S,  Slater_2S_3S},
    {  Slater_3S_1S,  Slater_3S_2S,  Slater_3S_3S}
};

t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]) = {
    {  DSlater_1S_1S,  DSlater_1S_2S,  DSlater_1S_3S},
    {  DSlater_2S_1S,  DSlater_2S_2S,  DSlater_2S_3S},
    {  DSlater_3S_1S,  DSlater_3S_2S,  DSlater_3S_3S}
};

t_slater_NS_func (*Slater_NS[SLATER_MAX]) = {
    Nuclear_1S,  Nuclear_2S,  Nuclear_3S
};

t_slater_NS_func (*DSlater_NS[SLATER_MAX]) = {
    DNuclear_1S,  DNuclear_2S,  DNuclear_3S
};

#endif
/* HAVE_LIBCLN */

extern "C" double Coulomb_SS(double r, int i, int j, double xi, double xj)
{
#ifdef HAVE_LIBCLN
    cl_R cr, cxi, cxj, cS;

    if ((i > SLATER_MAX) || (j > SLATER_MAX))
    {
        fprintf(stderr, "Slater-Slater integral %d %d not supported.\n", i, j);
        exit(1);
    }
    cxi = my_ftoa(xi);
    cxj = my_ftoa(xj);
    cr  = my_ftoa(r);
    if ((i > 0) && (j > 0))
    {
        cS = Slater_SS[i-1][j-1](cr, cxi, cxj);
    }
    else if (j > 0)
    {
        if (r == 0)
        {
            cS = cxj/j;
        }
        else
        {
            cS = Slater_NS[j-1](cr, cxj);
        }
    }
    else if (i > 0)
    {
        if (r == 0)
        {
            cS = cxi/i;
        }
        else
        {
            cS = Slater_NS[i-1](cr, cxi);
        }
    }
    else
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
    return double_approx(cS);
    
#else
    /* NOT HAVE_LIBCLN */
    
    double S = 0;
    if ((i > SLATER_MAX) || (j > SLATER_MAX))
    {
        fprintf(stderr, "Slater-Slater integral %d %d not supported.\n", i, j);
        exit(1);
    }
    if ((i > 0) && (j > 0))
    {
        S = Slater_SS[i-1][j-1](r, xi, xj);
    }
    else if (j > 0)
    {
        if (r == 0)
        {
            S = xj/j;
        }
        else
        {
            S = Slater_NS[j-1](r, xj);
        }
    }
    else if (i > 0)
    {
        if (r == 0)
        {
            S = xi/i;
        }
        else
        {
            S = Slater_NS[i-1](r, xi);
        }
    }
    else
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
    return S;
#endif
}

extern "C" double Nuclear_SS(double r, int i, double xi)
{
#ifdef HAVE_LIBCLN
    cl_R cr, cxi, cxj, cS;

    if (xi == 0)
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
    else if (r == 0)
    {
        return xi/i;
    }
    else if (i <= 0)
    {
        return 1/r;
    }
    else
    {
        cxi = my_ftoa(xi);
        cr  = my_ftoa(r);
        cS  = Slater_NS[i-1](cr, cxi);
        return double_approx(cS);
    }
#else
    double S = 0;
    if (xi == 0)
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
    else if (r == 0)
    {
        return xi/i;
    }
    else if (i <= 0)
    {
        return 1/r;
    }
    else
    {
        S  = Slater_NS[i-1](r, xi);
        return S;
    }
#endif
}

extern "C" double DCoulomb_SS(double r, int i, int j, double xi, double xj)
{
#ifdef HAVE_LIBCLN
    cl_R cr, cxi, cxj, cS;

    if ((i > SLATER_MAX) || (j > SLATER_MAX))
    {
        fprintf(stderr, "Slater-Slater integral %d %d not supported.\n", i, j);
        exit(1);
    }
    cxi = my_ftoa(xi);
    cxj = my_ftoa(xj);
    cr  = my_ftoa(r);
    if ((i > 0) && (j > 0))
    {
        cS = DSlater_SS[i-1][j-1](cr, cxi, cxj);
    }
    else if (j > 0)
    {
        if (r == 0)
        {
            cS = ZERO;
        }
        else
        {
            cS = DSlater_NS[j-1](cr, cxj);
        }
    }
    else if (i > 0)
    {
        if (r == 0)
        {
            cS = ZERO;
        }
        else
        {
            cS = DSlater_NS[i-1](cr, cxi);
        }
    }
    else
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1/(r*r);
        }
    }
    return double_approx(cS);
#else
    double S = 0;
    
    if ((i > SLATER_MAX) || (j > SLATER_MAX))
    {
        fprintf(stderr, "Slater-Slater integral %d %d not supported.\n", i, j);
        exit(1);
    }
    if ((i > 0) && (j > 0))
    {
        S = DSlater_SS[i-1][j-1](r, xi, xj);
    }
    else if (j > 0)
    {
        if (r == 0)
        {
            S = 0;
        }
        else
        {
            S = DSlater_NS[j-1](r, xj);
        }
    }
    else if (i > 0)
    {
        if (r == 0)
        {
            S = 0;
        }
        else
        {
            S = DSlater_NS[i-1](r, xi);
        }
    }
    else
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1/(r*r);
        }
    }
    return S;
#endif
}

extern "C" double DNuclear_SS(double r, int i, double xi)
{
#ifdef HAVE_LIBCLN
    cl_R cr, cxi, cxj, cS;

    if (i > SLATER_MAX)
    {
        fprintf(stderr, "Slater-Nuclear integral %d not supported.\n", i);
        exit(1);
    }
    if (r == 0)
    {
        return 0;
    }
    else
    {
        if ((xi == 0) || (i <= 0))
        {
            return -1/(r*r);
        }
        else
        {
            cxi = my_ftoa(xi);
            cr  = my_ftoa(r);
            cS  = DSlater_NS[i-1](cr, cxi);
            return double_approx(cS);
        }
    }
#else
    double S = 0;

    if (i > SLATER_MAX)
    {
        fprintf(stderr, "Slater-Nuclear integral %d not supported.\n", i);
        exit(1);
    }
    if (r == 0)
    {
        return 0;
    }
    else
    {
        if ((xi == 0) || (i <= 0))
        {
            return -1/(r*r);
        }
        else
        {
            S  = DSlater_NS[i-1](r, xi);
            return S;
        }
    }
#endif
}

#ifdef HAVE_LIBCLN
cl_R Power(cl_R a, int b)
{
    int  minus = 0;
    cl_R cP;

    if (b < 0)
    {
        minus = 1;
        b     = -b;
    }
    if (b == 0)
    {
        return ONE;
    }
    if (a == ZERO)
    {
        return ZERO;
    }
    if ((b % 2) == 0)
    {
        cP = Power(a*a, b/2);
        if (minus)
        {
            return ONE/cP;
        }
        else
        {
            return cP;
        }
    }
    else if ((b % 2) == 1)
    {
        cP = a*Power(a*a, b/2);
        if (minus)
        {
            return ONE/cP;
        }
        else
        {
            return cP;
        }
    }
    return ZERO;
}
#endif
