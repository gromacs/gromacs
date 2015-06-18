/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal\file
 * \brief Implements utilities for initializing long range correction tables
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "longrangetables.h"

#include <cmath>

#include <algorithm>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/tables/splinetable.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

void v_q_ewald_lr(const std::vector<double> &params,
                  double r, double *e, double *f)
{
    double br = params[0] * r;
    if (r == 0)
    {
        *e = params[0]*2/sqrt(M_PI);
        *f = 0;
    }
    else
    {
        *e = std::erf(br)/r;
        /* Force derived using Mathematica from
         * -D[Erf[params r]/r, r] (which is 0 in the limit r -> 0)
         */
        *f = std::erf(br)/(r*r) - 2*params[0]*exp(-br*br)/(r*sqrt(M_PI));
    }
}

/*! \brief
 * Compute the long range energy and force for Ewald Lennard Jones (r^{-6} only)
 *
 * This function is completely analytical which may lead to instabilities
 * in the returned energy and force. See the next fucntion for more details.
 * \param[in]  params Parameters to the long range function
 * \param[in]  r    Distance
 * \param[out] e    Energy
 * \param[out] f    Force
 */
static void v_lj_ewald_lr_analytical(const std::vector<double> &params,
                                     double r, double *e, double *f)
{
    if (r == 0)
    {
        *e = gmx::power6(params[0])/6.0;
        *f = 0;
    }
    else
    {
        double br      = params[0]*r;
        double br2     = br*br;
        double br4     = br2*br2;
        double br6     = br4*br2;
        double r6      = gmx::power6(r);
        double expmbr2 = std::exp(-br2);
        double gmxem   = std::expm1(-br2);
        /* In order to preserve (obtain) accuracy it is necessary to use
         * the expm1 function to replace (exp(x)-1) which becomes close to 0
         * as x -> 0.
         * Original energy expression:
         * *e  = (1.0 - exp(-br2)*(1 + br2 + 0.5*br4))/r6;
         */
        *e            = -(gmxem + expmbr2*(br2 + 0.5*br4))/r6;
        /* Force derived using Mathematica from
         * Simplify[-D[(1 - Exp[-b^2 r^2]*(1 + b^2 r^2 + (b^4 r^4)/2))/(r^6), r]]
         * See above about expm1, original expression:
         * *f  = -(expmbr2*(6 - 6/expmbr2 + 6*br2 + 3*br4 + br2*br4))/(r6*r);
         */
        *f            = -(6*(gmxem + br2*expmbr2) + (3*br4 + br6)*expmbr2)/(r6*r);
    }
}

void v_lj_ewald_lr(const std::vector<double> &params,
                   double r, double *e, double *f)
{
    v_lj_ewald_lr_analytical(params, r, e, f);
    double rmin = 0.1, rmax = 0.2;
    if ((r == 0) || (r > rmax))
    {
        return;
    }

    /* Start from Mathematica code:
     * LRcorrDispV[r_] := (1 - Exp[-b^2 r^2]*(1 + b^2 r^2 + (b^4 r^4)/2))/(r^6)
     * ApproxV[r_] = Normal[Series[LRcorrDispV[r], {r, 0, 16}]]
     * Substituting x for b^2r^2 this gives
     * b^6(1/6 - x/8 + x^2/20 - x^3/72 + x^4/336 - x^5/1920 + x^6/12960 - x^7/100800 + x^8/887040)
     *
     * LRcorrDispF[r_] = -Simplify[D[LRcorrDispV[r], r]]
     * ApproxF[r_] = Normal[Series[LRcorrDispF[r], {r, 0, 16}]]
     * using the same substitution, this gives
     * b^8r(1/4 - x/5 + x^2/12 - x^3/42 + x^4/192 - x^5/1080 + x^6/7200 - x^7/55440)
     */
    double x  = pow(params[0]*r, 2);
    double ea = pow(params[0], 6)*(pow(x, 7)*(x/887040.0 - 1/100800.0) + pow(x, 5)*(x/12960.0 - 1/1920.0) +
                                   pow(x, 3)*(x/336.0 - 1/72.0) + x*(x/20.0 - 1/8.0) + 1/6.0);
    double fa = pow(params[0], 8)*r*(pow(x, 6)*(-x/55440.0 + 1/7200.0) + pow(x, 4)*(-x/1080.0 + 1/192.0) +
                                     pow(x, 2)*(-x/42.0 + 1/12.0) - x/5.0 + 1/4.0);
    /* Now do linear interpolation between approximation and
     * analytical function.
     */
    if (r <= rmin)
    {
        *e = ea;
        *f = fa;
    }
    else
    {
        double dr = (r-rmin)/(rmax-rmin);
        *e = (1-dr)*ea + dr* *e;
        *f = (1-dr)*fa + dr* *f;
    }
}

/*! \brief Return spacing for tables.
 *
 * Returns the spacing for a function using the maximum of
 * the third derivative, x_scale (unit 1/length)
 * and function tolerance.
 * \param[in] third_deriv_max A number
 * \param[in] x_scale         The distance between table points
 * \param[in] func_tol        Maximum precision
 * \return Table spacing
 */
static double spline3_table_scale(double third_deriv_max,
                                  double x_scale,
                                  double func_tol)
{
    double deriv_tol;
    double sc_deriv, sc_func;

    /* Force tolerance: single precision accuracy */
    deriv_tol = GMX_FLOAT_EPS;
    sc_deriv  = sqrt(third_deriv_max/(6*4*deriv_tol*x_scale))*x_scale;

    /* Don't try to be more accurate on energy than the precision */
    func_tol  = std::max(func_tol, static_cast<double>(GMX_REAL_EPS));
    sc_func   = std::pow(third_deriv_max/(6*12*sqrt(3.0)*func_tol), 1.0/3.0)*x_scale;

    return std::max(sc_deriv, sc_func);
}

double longRangeCorrectionTableScale(int    type_coul,
                                     double ewaldcoeff_coul,
                                     double r_coul,
                                     int    type_vdw,
                                     double ewaldcoeff_vdw,
                                     double r_vdw)
{
    double sc;

    sc = 0;

    if (type_coul == eelEWALD || EEL_PME(type_coul))
    {
        double erf_x_d3 = 1.0522; /* max of (erf(x)/x)''' */
        double etol;
        double sc_q;

        /* Energy tolerance: 0.1 times the cut-off jump */
        etol  = 0.1*std::erfc(ewaldcoeff_coul*r_coul);

        sc_q  = spline3_table_scale(erf_x_d3, ewaldcoeff_coul, etol);

        if (debug)
        {
            fprintf(debug, "Ewald Coulomb quadratic spline table spacing: %f 1/nm\n", 1/sc_q);
        }

        sc    = std::max(sc, sc_q);
    }

    if (EVDW_PME(type_vdw))
    {
        double func_d3 = 0.42888; /* max of (x^-6 (1 - exps(-x^2)(1+x^2+x^4/2)))''' */
        double xrc2, etol;
        double sc_lj;

        /* Energy tolerance: 0.1 times the cut-off jump */
        xrc2  = gmx::square(ewaldcoeff_vdw*r_vdw);
        etol  = 0.1*exp(-xrc2)*(1 + xrc2 + xrc2*xrc2/2.0);

        sc_lj = spline3_table_scale(func_d3, ewaldcoeff_vdw, etol);

        if (debug)
        {
            fprintf(debug, "Ewald LJ quadratic spline table spacing: %f 1/nm\n", 1/sc_lj);
        }

        sc = std::max(sc, sc_lj);
    }

    return sc;
}

unsigned int longRangeCorrectionTableSize(double rtab,
                                          double scale)
{
    return static_cast<unsigned int>(rtab*scale)+1;
}
