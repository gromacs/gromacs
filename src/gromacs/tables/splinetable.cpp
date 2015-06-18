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
/*! \brief
 * Implements class for spline interpolation tables
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 */
#include "gmxpre.h"

#include "splinetable.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/snprintf.h"

/*! \brief Check that the force is the negative derivative of the energy
 *
 * \param[in] ntable    The number of data points
 * \param[in] dx        The distance between points
 * \param[in] v         The energy array
 * \param[in] f         The force array
 * \param[in] rel_toler Relative tolerance
 * \throws If the average relative error is larger than the tolerance
 */
static void check_force(unsigned int  ntable,
                        double        dx,
                        const double *v,
                        const double *f,
                        double        rel_toler)
{
    if (ntable < 4)
    {
        return;
    }
    double       fsum_err = 0;
    unsigned int nsum     = 0;
    for (unsigned int i = 1; (i < ntable-1); i++)
    {
        double   dv = (v[i+1]-v[i-1])/(2*dx);
        double   df = (f[i]+dv);
        if (f[i] != 0)
        {
            fsum_err   += fabs(df/f[i]);
            nsum++;
        }
    }
    double faver = fsum_err/nsum;
    char   buf[256];
    snprintf(buf, sizeof(buf), "Force input to tables deviates too much from derivative of energy.\nN = %u, average relative error = %g tolerance = %g.", ntable, faver, rel_toler);
    GMX_RELEASE_ASSERT((faver <= rel_toler), buf);
}

namespace gmx
{

QuadraticSplineTable::QuadraticSplineTable(unsigned int  ntable,
                                           double        dx,
                                           const double *table_in_v,
                                           const double *table_in_f)
{
    unsigned int stride;

    GMX_RELEASE_ASSERT((table_in_v != NULL), "No table input passed");
    GMX_RELEASE_ASSERT((ntable >= 2),
                       "Can not make a spline table with less than 2 points");

    if (NULL != table_in_f)
    {
        stride = 1;
        check_force(ntable, dx, table_in_v, table_in_f, 0.001);
    }
    else
    {
        stride    = 2;
        ntable   /= 2;
        dx       *= 2;
    }
    tabSpacing_ = dx;
    tabSize_    = ntable;

    fillTable(stride, table_in_v, table_in_f);
}

QuadraticSplineTable::QuadraticSplineTable(unsigned int                    ntable,
                                           double                          dx,
                                           interaction_potential_function  v_ana,
                                           const double                   *params)
{
    std::vector<double> table_v;
    std::vector<double> table_f;

    GMX_RELEASE_ASSERT((v_ana != NULL), "No analytical function passed");

    GMX_RELEASE_ASSERT((ntable >= 2),
                       "Can not make a spline table with less than 2 points");
    tabSpacing_ = dx;
    tabSize_    = ntable;
    for (unsigned int i = 0; (i < tabSize_); i++)
    {
        double x, e, f;
        x = i*tabSpacing_;
        v_ana(params, x, &e, &f);
        table_v.push_back(e);
        table_f.push_back(f);
    }
    fillTable(1, &(table_v[0]), &(table_f[0]));
}

void QuadraticSplineTable::fillTable(unsigned int  stride,
                                     const double *table_in_v,
                                     const double *table_in_f)
{
    int     i_inrange;
    double  tab_max;
    double  dc, dc_new;
    bool    bOutOfRange;
    double  v_r0, v_r1, v_inrange, vi, a0, a1, a2dx;
    /* Allocate memory */
    tableF_.resize(tabSize_);
    tableV_.resize(tabSize_);
    tableFDV0_.resize(4*tabSize_);

    /* We need some margin to be able to divide table values by r
     * in the kernel and also to do the integration arithmetics
     * without going out of range. Furthemore, we divide by tabSpacing_ below.
     */
    tab_max = GMX_REAL_MAX*0.0001;

    /* This function produces a table with:
     * maximum energy error: V'''/(6*12*sqrt(3))*tabSpacing_^3
     * maximum force error:  V'''/(6*4)*tabSpacing_^2
     * The rms force error is the max error times 1/sqrt(5)=0.45.
     */
    bOutOfRange = false;
    i_inrange   = tabSize_;
    v_inrange   = 0;
    dc          = 0;
    for (unsigned int i = tabSize_-1; true; i--)
    {
        v_r0 = table_in_v[i*stride];

        if (!bOutOfRange)
        {
            i_inrange = i;
            v_inrange = v_r0;

            vi = v_r0;
        }
        else
        {
            /* Linear continuation for the last point in range */
            vi = v_inrange - dc*(i - i_inrange)*tabSpacing_;
        }

        tableV_[i] = vi;

        if (i == 0)
        {
            break;
        }
        /* Get the potential at table point i-1 */
        v_r1 = table_in_v[(i - 1)*stride];

        if (v_r1 != v_r1 || v_r1 < -tab_max || v_r1 > tab_max)
        {
            bOutOfRange = true;
        }

        if (!bOutOfRange)
        {
            /* Calculate the average second derivative times tabSpacing_ over interval i-1 to i */
            if (table_in_f == NULL)
            {
                /* Use the function values at the end points and the middle */
                double v_mid;

                v_mid = table_in_v[i*stride - stride/2];
                a2dx  = (v_r0 + v_r1 - 2*v_mid)/(0.25*tabSpacing_);
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
            dc = (v_r0 - v_r1)/tabSpacing_ + 0.5*a2dx;
        }

        if (i == tabSize_ - 1)
        {
            /* Fill the table with the force, minus the derivative of the spline */
            tableF_[i] = -dc;
        }
        else
        {
            /* tab[i] will contain the average of the splines over the two intervals */
            tableF_[i] += -0.5*dc;
        }

        if (!bOutOfRange)
        {
            /* Make spline s(x) = a0 + a1*(x - xr) + 0.5*a2*(x - xr)^2
             * matching the potential at the two end points
             * and the derivative dc at the end point xr.
             */
            a0   = v_r0;
            a1   = dc;
            a2dx = (a1*tabSpacing_ + v_r1 - a0)*2/tabSpacing_;

            /* Set dc to the derivative at the next point */
            dc_new = a1 - a2dx;

            if (dc_new != dc_new || dc_new < -tab_max || dc_new > tab_max)
            {
                bOutOfRange = true;
            }
            else
            {
                dc = dc_new;
            }
        }

        tableF_[(i-1)] = -0.5*dc;
    }
    /* Currently the last value only contains half the force: double it */
    tableF_[0] *= 2;

    /* Copy to FDV0 table too. */
    for (unsigned int i = 0; i < tabSize_-1; i++)
    {
        tableFDV0_[4*i]     = tableF_[i];
        tableFDV0_[4*i+1]   = tableF_[i+1]-tableF_[i];
        tableFDV0_[4*i+2]   = tableV_[i];
        tableFDV0_[4*i+3]   = 0.0;
    }
    tableFDV0_[4*(tabSize_-1)]    = tableF_[(tabSize_-1)];
    // We are missing the next force hence this is incorrect.
    tableFDV0_[4*(tabSize_-1)+1]  = -tableF_[(tabSize_-1)];
    tableFDV0_[4*(tabSize_-1)+2]  = tableV_[(tabSize_-1)];
    tableFDV0_[4*(tabSize_-1)+3]  = 0.0;
    // To fix the problem we decrease the table size here.
    tabSize_ -= 1;
}

double QuadraticSplineTable::force(double r) const
{
    double       dx   = spacing();
    GMX_RELEASE_ASSERT((r >= 0), "Distance r should be >= 0");
    unsigned int itab = static_cast<int>(std::floor(r/dx));
    GMX_RELEASE_ASSERT((itab < size()), "Distance r too large for table");
    double       dr   = r/dx - itab;
    double       f    = tableFDV0_[4*itab]+dr*tableFDV0_[4*itab+1];
    return f;
}

double QuadraticSplineTable::energy(double r) const
{
    double       dx   = spacing();
    GMX_RELEASE_ASSERT((r >= 0), "Distance r should be >= 0");
    unsigned int itab = static_cast<int>(std::floor(r/dx));
    GMX_RELEASE_ASSERT((itab < size()), "Distance r too large for table");
    double       dr   = r - dx*itab;
    double       v    = (tableFDV0_[4*itab+2] - (dr*(tableFDV0_[4*itab] + (dr/(2*dx))*tableFDV0_[4*itab+1])));
    return v;
}

} // namespace gmx
