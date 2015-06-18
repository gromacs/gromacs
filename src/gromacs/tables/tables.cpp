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
#include "gmxpre.h"

#include "tables.h"

#include <cmath>

#include <algorithm>

#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

double v_q_ewald_lr(const double *beta, double r)
{
    if (r == 0)
    {
        return (*beta)*2/sqrt(M_PI);
    }
    else
    {
        return gmx_erfd((*beta)*r)/r;
    }
}

double v_lj_ewald_lr(const double *beta, double r)
{
    double br, br2, br4, r6, factor;
    if (r == 0)
    {
        return std::pow((*beta), 6.0)/6.0;
    }
    else
    {
        br     = (*beta)*r;
        br2    = br*br;
        br4    = br2*br2;
        r6     = std::pow(r, 6.0);
        factor = (1.0 - exp(-br2)*(1 + br2 + 0.5*br4))/r6;

        return factor;
    }
}

static void table_spline3_fill(std::vector<double>            &table_f,
                               std::vector<double>            &table_v,
                               std::vector<double>            &table_fdv0,
                               unsigned int                    ntab,
                               double                          dx,
                               interaction_potential_function  v_ana,
                               const double                   *params,
                               const double                   *table_in_v,
                               const double                   *table_in_f,
                               unsigned int                    table_in_size,
                               unsigned int                    stride)
{
    double       tab_max;
    int          i_inrange;
    double       dc, dc_new;
    gmx_bool     bOutOfRange;
    double       v_r0, v_r1, v_inrange, vi, a0, a1, a2dx;
    double       x_r0;

    GMX_RELEASE_ASSERT(!(v_ana != NULL && table_in_v != NULL),
                       "Both analytical function and table input passed");

    if (table_in_v != NULL)
    {
        char buf[256];
        GMX_RELEASE_ASSERT((stride > 0),
                           "The stride for the table input should be positive");

        snprintf(buf, sizeof(buf),
                 "For filling a quadratic spline table with tabulated potential input only (i.e. no force input), the stride has to be a mulitple of 2, not odd (%u)", stride);
        GMX_RELEASE_ASSERT(!(table_in_f == NULL && stride % 2 != 0), buf);

        snprintf(buf, sizeof(buf),
                 "The table input size (%u) is smaller than the requested output table size (%u) times the stride (%u)",
                 table_in_size, ntab, stride);
        GMX_RELEASE_ASSERT((table_in_size >= ntab*stride), buf);
    }

    GMX_RELEASE_ASSERT((ntab >= 2),
                       "Can not make a spline table with less than 2 points");

    /* Allocate memory */
    table_f.resize(ntab);
    table_v.resize(ntab);
    table_fdv0.resize(4*ntab);

    /* We need some margin to be able to divide table values by r
     * in the kernel and also to do the integration arithmetics
     * without going out of range. Furthemore, we divide by dx below.
     */
    tab_max = GMX_REAL_MAX*0.0001;

    /* This function produces a table with:
     * maximum energy error: V'''/(6*12*sqrt(3))*dx^3
     * maximum force error:  V'''/(6*4)*dx^2
     * The rms force error is the max error times 1/sqrt(5)=0.45.
     */

    bOutOfRange = FALSE;
    i_inrange   = ntab;
    v_inrange   = 0;
    dc          = 0;
    for (unsigned int i = ntab-1; true; i--)
    {
        x_r0 = i*dx;

        if (v_ana != NULL)
        {
            v_r0 = (*v_ana)(params, x_r0);
        }
        else
        {
            v_r0 = table_in_v[i*stride];
        }

        if (!bOutOfRange)
        {
            i_inrange = i;
            v_inrange = v_r0;

            vi = v_r0;
        }
        else
        {
            /* Linear continuation for the last point in range */
            vi = v_inrange - dc*(i - i_inrange)*dx;
        }

        table_v[i] = vi;

        if (i == 0)
        {
            break;
        }
        if (i < 5)
        {
            printf("i = %u\n", i);
        }
        /* Get the potential at table point i-1 */
        if (v_ana != NULL)
        {
            v_r1 = (*v_ana)(params, (i-1)*dx);
        }
        else
        {
            v_r1 = table_in_v[(i - 1)*stride];
        }

        if (v_r1 != v_r1 || v_r1 < -tab_max || v_r1 > tab_max)
        {
            bOutOfRange = TRUE;
        }

        if (!bOutOfRange)
        {
            /* Calculate the average second derivative times dx over interval i-1 to i */
            if (v_ana != NULL || table_in_f == NULL)
            {
                /* Use the function values at the end points and the middle */
                double v_mid;

                if (v_ana != NULL)
                {
                    v_mid = (*v_ana)(params, x_r0 - 0.5*dx);
                }
                else
                {
                    v_mid = table_in_v[i*stride - stride/2];
                }
                a2dx = (v_r0 + v_r1 - 2*v_mid)/(0.25*dx);
            }
            else
            {
                /* Use the forces at the end points */
                a2dx = -table_in_f[i*stride] + table_in_f[(i - 1)*stride];
            }
            /* Set the derivative of the spline to match the difference in potential
             * over the interval plus the average effect of the quadratic term.
             * This is the essential step for minimizing the error in the force.
             */
            dc = (v_r0 - v_r1)/dx + 0.5*a2dx;
        }

        if (i == ntab - 1)
        {
            /* Fill the table with the force, minus the derivative of the spline */
            table_f[i] = -dc;
        }
        else
        {
            /* tab[i] will contain the average of the splines over the two intervals */
            table_f[i] += -0.5*dc;
        }

        if (!bOutOfRange)
        {
            /* Make spline s(x) = a0 + a1*(x - xr) + 0.5*a2*(x - xr)^2
             * matching the potential at the two end points
             * and the derivative dc at the end point xr.
             */
            a0   = v_r0;
            a1   = dc;
            a2dx = (a1*dx + v_r1 - a0)*2/dx;

            /* Set dc to the derivative at the next point */
            dc_new = a1 - a2dx;

            if (dc_new != dc_new || dc_new < -tab_max || dc_new > tab_max)
            {
                bOutOfRange = TRUE;
            }
            else
            {
                dc = dc_new;
            }
        }

        table_f[(i-1)] = -0.5*dc;
    }
    /* Currently the last value only contains half the force: double it */
    table_f[0] *= 2;

    /* Copy to FDV0 table too. */
    for (unsigned int i = 0; i < ntab-1; i++)
    {
        table_fdv0[4*i]     = table_f[i];
        table_fdv0[4*i+1]   = table_f[i+1]-table_f[i];
        table_fdv0[4*i+2]   = table_v[i];
        table_fdv0[4*i+3]   = 0.0;
    }
    table_fdv0[4*(ntab-1)]    = table_f[(ntab-1)];
    table_fdv0[4*(ntab-1)+1]  = -table_f[(ntab-1)];
    table_fdv0[4*(ntab-1)+2]  = table_v[(ntab-1)];
    table_fdv0[4*(ntab-1)+3]  = 0.0;
}

namespace gmx
{
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

void InteractionTables::initTabScaleLength(int    type_coul,
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
        etol  = 0.1*gmx_erfc(ewaldcoeff_coul*r_coul);

        sc_q  = spline3_table_scale(erf_x_d3, ewaldcoeff_coul, etol);

        if (debug)
        {
            fprintf(debug, "Ewald Coulomb quadratic spline table spacing: %f 1/nm\n", 1/sc_q);
        }

        sc    = std::max(sc, sc_q);
    }

    if (EVDW_PME(type_vdw))
    {
        double func_d3 = 0.42888; /* max of (x^-6 (1 - exp(-x^2)(1+x^2+x^4/2)))''' */
        double xrc2, etol;
        double sc_lj;

        /* Energy tolerance: 0.1 times the cut-off jump */
        xrc2  = sqr(ewaldcoeff_vdw*r_vdw);
        etol  = 0.1*exp(-xrc2)*(1 + xrc2 + xrc2*xrc2/2.0);

        sc_lj = spline3_table_scale(func_d3, ewaldcoeff_vdw, etol);

        if (debug)
        {
            fprintf(debug, "Ewald LJ quadratic spline table spacing: %f 1/nm\n", 1/sc_lj);
        }

        sc = std::max(sc, sc_lj);
    }

    tabScale_  = 1.0/sc;
    tabLength_ = std::max(r_coul, r_vdw);
    tabSize_   = 1+static_cast<int>(tabLength_/tabScale_);
}

InteractionTables::InteractionTables(int                             type_coul,
                                     double                          ewaldcoeff_coul,
                                     double                          r_coul,
                                     double                         *params_coul,
                                     interaction_potential_function  v_coul,
                                     const double                   *table_v_coul,
                                     const double                   *table_f_coul,
                                     int                             type_vdw,
                                     double                          ewaldcoeff_vdw,
                                     double                          r_vdw,
                                     double                         *params_vdw,
                                     interaction_potential_function  v_vdw,
                                     const double                   *table_v_vdw,
                                     const double                   *table_f_vdw,
                                     unsigned int                    table_in_length,
                                     unsigned int                    table_in_stride)
{
    initTabScaleLength(type_coul, ewaldcoeff_coul, r_coul,
                       type_vdw,  ewaldcoeff_vdw,  r_vdw);

    if (type_coul == eelEWALD || EEL_PME(type_coul))
    {
        table_spline3_fill(tabCoulF_, tabCoulV_, tabCoulFDV0_,
                           tabSize_, tabScale_, v_coul, params_coul,
                           table_v_coul, table_f_coul,
                           table_in_length, table_in_stride);
    }
    if (EVDW_PME(type_vdw))
    {
        table_spline3_fill(tabVdwF_, tabVdwV_, tabVdwFDV0_,
                           tabSize_, tabScale_, v_vdw, params_vdw,
                           table_v_vdw, table_f_vdw,
                           table_in_length, table_in_stride);
    }
}

} // namespace gmx
