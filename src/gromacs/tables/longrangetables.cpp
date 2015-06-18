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

#include "longrangetables.h"

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
#include "gromacs/tables/splinetable.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

void v_q_ewald_lr(const double *beta, double r, double *e, double *f)
{
    double br = *beta * r;
    if (r == 0)
    {
        *e = (*beta)*2/sqrt(M_PI);
        *f = 0;
    }
    else
    {
        *e = gmx_erfd(br)/r;
        /* Force derived using Mathematica from
         * -D[Erf[beta r]/r, r] (which is 0 in the limit r -> 0)
         */
        *f = gmx_erfd(br)/(r*r) - 2*(*beta)*exp(-br*br)/(r*sqrt(M_PI));
    }
}

void v_lj_ewald_lr(const double *beta, double r, double *e, double *f)
{
    if (r == 0)
    {
        *e = pow((*beta), 6.0)/6.0;
        *f = 0;
    }
    else
    {
        double br     = (*beta)*r;
        double br2    = br*br;
        double br4    = br2*br2;
        double r6     = pow(r, 6.0);
        double expbr2 = exp(-br2);
        /* In order to preserve (obtain) accuracy it is necessary to use
         * the gmx_expm1 function to replace (1-exp(x)) which becomes close to 0
         * as x -> 0.
         * Original energy expression:
         * *e  = (1.0 - expbr2*(1 + br2 + 0.5*br4))/r6;
         */
        *e            = -(gmx_expm1(-br2) + expbr2*(br2 + 0.5*br4))/r6;
        /* Force derived using Mathematica from
         * Simplify[-D[(1 - Exp[-b^2 r^2]*(1 + b^2 r^2 + (b^4 r^4)/2))/(r^6), r]]
         * See above about gmx_expm1, original expression:
         * *f  = -(expbr2*(6 - 6/expbr2 + 6*br2 + 3*br4 + br2*br4))/(r6*r);
         */
        *f            = -(6*gmx_expm1(-br2) + (6*br2 + (3*br4 + br2*br4))*expbr2)/(r6*r);
    }
}

static gmx::SplineTable *table_spline3_fill(unsigned int                    ntab,
                                            double                          dx,
                                            interaction_potential_function  v_ana,
                                            const double                   *params)
{
    std::vector<double> table_v;
    std::vector<double> table_f;

    GMX_RELEASE_ASSERT((v_ana != NULL), "No analytical function passed");

    GMX_RELEASE_ASSERT((ntab >= 2),
                       "Can not make a spline table with less than 2 points");
    for (unsigned int i = 0; (i < ntab); i++)
    {
        double x, e, f;
        x = i*dx;
        v_ana(params, x, &e, &f);
        table_v.push_back(e);
        table_f.push_back(f);
    }
    return new gmx::SplineTable(ntab, dx, &(table_v[0]), &(table_f[0]));
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

double LongRangeCorrectionTables::computeSpacing(int    type_coul,
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
        double func_d3 = 0.42888; /* max of (x^-6 (1 - exps(-x^2)(1+x^2+x^4/2)))''' */
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

    return 1.0/sc;
}

double LongRangeCorrectionTables::spacing() const
{
    if (NULL != tableCoulomb_)
    {
        return tableCoulomb_->spacing();
    }
    else if (NULL != tableVdw_)
    {
        return tableVdw_->spacing();
    }
    return 0.0;
}

unsigned int LongRangeCorrectionTables::size()
{
    if (NULL != tableCoulomb_)
    {
        return tableCoulomb_->size();
    }
    else if (NULL != tableVdw_)
    {
        return tableVdw_->size();
    }
    return 0.0;
}

LongRangeCorrectionTables::LongRangeCorrectionTables(int                             type_coul,
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
                                                     unsigned int                    table_length,
                                                     double                          table_dx,
                                                     double                          r_tab)
{
    double   my_dx     = computeSpacing(type_coul, ewaldcoeff_coul, r_coul,
                                        type_vdw,  ewaldcoeff_vdw,  r_vdw);
    double   tabLength = r_tab + std::max(r_coul, r_vdw);
    unsigned tabSize   = 1+static_cast<int>(tabLength/my_dx);

    tableCoulomb_ = NULL;
    tableVdw_     = NULL;
    if (type_coul == eelEWALD || EEL_PME(type_coul))
    {
        if (NULL != table_v_coul)
        {
            tableCoulomb_ = new SplineTable(table_length, table_dx,
                                            table_v_coul, table_f_coul);
        }
        else
        {
            tableCoulomb_ = table_spline3_fill(tabSize, my_dx, v_coul, params_coul);
        }
    }
    if (EVDW_PME(type_vdw))
    {
        if (NULL != table_v_vdw)
        {
            tableVdw_   = new SplineTable(table_length, table_dx,
                                          table_v_vdw, table_f_vdw);
        }
        else
        {
            tableVdw_ = table_spline3_fill(tabSize, my_dx, v_vdw, params_vdw);
        }
    }
    if ((NULL != tableCoulomb_) && (NULL != tableVdw_))
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "Spacing for Coulomb %g and Vdw %g\n",
                 tableCoulomb_->spacing(),
                 tableVdw_->spacing());
        GMX_RELEASE_ASSERT((tableVdw_->spacing() == tableCoulomb_->spacing()),
                           buf);
    }
}

LongRangeCorrectionTables::~LongRangeCorrectionTables()
{
    if (NULL != tableCoulomb_)
    {
        delete tableCoulomb_;
    }
    if (NULL != tableVdw_)
    {
        delete tableVdw_;
    }
}

} // namespace gmx
