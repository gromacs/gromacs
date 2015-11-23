/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares data structure and utilities for electric fields
 *
 * \inpublicapi
 * \ingroup module_mdtypes
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "electricfield.h"

#include <cmath>

#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

void ElectricField::setFieldTerm(int dim, real a, real omega, real t0, real sigma)
{
    range_check(dim, 0, DIM);
    efield_[dim].setField(a, omega, t0, sigma);
}

real ElectricField::field(int dim, real t) const
{
    range_check(dim, 0, DIM);
    if (efield_[dim].sigma() > 0)
    {
        real t0 = efield_[dim].t0();
        return efield_[dim].a() * (cos(efield_[dim].omega()*(t-t0))*
                                   exp(-sqr(t-t0)/(2.0*sqr(efield_[dim].sigma()))));
    }
    else
    {
        return efield_[dim].a() * cos(efield_[dim].omega()*t);
    }
}

bool ElectricField::applyField() const
{
    return (efield_[XX].a() != 0 ||
            efield_[YY].a() != 0 ||
            efield_[ZZ].a() != 0);
}

void ElectricField::printComponents(FILE *fp, double t)
{
    try
    {
        fprintf(fp, "%10g  %10g  %10g  %10g #FIELD\n", t,
                field(XX, t), field(YY, t), field(ZZ, t));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void calc_f_el(int  start, int homenr,
               real charge[], rvec f[],
               const ElectricField *efield,
               double t)
{
    for (int m = 0; (m < DIM); m++)
    {
        real Ext = FIELDFAC*efield->field(m, t);

        if (Ext != 0)
        {
            for (int i = start; (i < start+homenr); i++)
            {
                f[i][m] += charge[i]*Ext;
            }
        }
    }
}
