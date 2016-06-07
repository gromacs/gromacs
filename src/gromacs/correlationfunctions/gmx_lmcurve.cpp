/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal
 * \file
 * \brief
 * Defines a driver routine for lmfit, and a callback for it to use.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gmx_lmcurve.h"

#include <lmmin.h>
#include <lmstruct.h>

typedef struct {
    const double* t;
    const double* y;
    const double* dy;
    double (*f)(const double t, const double* par);
} lmcurve_data_struct;

//! Callback function used by lmmin
static void lmcurve_evaluate(
        const double* par, const int m_dat, const void* data, double* fvec,
        int* info)
{
    lmcurve_data_struct* D = (lmcurve_data_struct*)data;
    int                  i;
    for (i = 0; i < m_dat; i++)
    {
        double dy = D->dy[i];
        if (dy == 0)
        {
            dy = 1;
        }
        fvec[i] = (D->y[i] - D->f(D->t[i], par))/dy;
    }
    *info = 0;
}

void gmx_lmcurve(
        const int n_par, double* par, const int m_dat,
        const double* t, const double* y, const double *dy,
        double (*f)(double t, const double* par),
        const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_data_struct data = { t, y, dy, f };

    lmmin(n_par, par, m_dat, (const void*)&data, lmcurve_evaluate,
          control, status);
}
