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
 * \inlibraryapi
 * \ingroup module_applied_forces
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "electricfield.h"

#include <cmath>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/fileio/txtdump.h"
#include "gromacs/gmxlib/compare.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/warninp.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

void ElectricField::doTpxIO(t_fileio *fio, bool bRead)
{
    for (int j = 0; (j < DIM); j++)
    {
        int n = 1, nt = 1;
        gmx_fio_do_int(fio, n);
        gmx_fio_do_int(fio, nt);
        std::vector<real> aa, phi, at, phit;
        if (!bRead)
        {
            aa.push_back(a(j));
            phi.push_back(t0(j));
            at.push_back(omega(j));
            phit.push_back(sigma(j));
        }
        else
        {
            aa.resize(n+1);
            phi.resize(n+1);
            at.resize(nt+1);
            phit.resize(nt+1);
        }
        gmx_fio_ndo_real(fio, aa.data(),  n);
        gmx_fio_ndo_real(fio, phi.data(), n);
        gmx_fio_ndo_real(fio, at.data(),  nt);
        gmx_fio_ndo_real(fio, phit.data(), nt);
        if (bRead && n > 0 && nt > 0)
        {
            setFieldTerm(j, aa[0], phi[0], at[0], phit[0]);
        }
    }

}

void ElectricField::decodeMdp(int         dim,
                              const char *staticField,
                              const char *dynamicField,
                              warninp_t   wi)
{
    std::vector<std::string> sx  = gmx::splitString(staticField);
    if (sx.size() >= 3)
    {
        real aa, omega, t0, sigma;
        int  n;
        // Read old style mdp files
        try
        {
            char *pos;
            n  = std::strtol(sx[0].c_str(), &pos, 10);
            if (n > 1)
            {
                char warn_buf[STRLEN];
                sprintf(warn_buf, "Only one electric field term supported for each dimension");
                warning_error(wi, warn_buf);
            }
            aa = std::strtod(sx[1].c_str(), &pos);
            t0 = std::strtod(sx[2].c_str(), &pos);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

        std::vector<std::string> sxt = gmx::splitString(dynamicField);
        if (sxt.size() >= 3)
        {
            try
            {
                char *pos;
                n  = std::strtol(sx[0].c_str(), &pos, 10);
                if (n > 1)
                {
                    char warn_buf[STRLEN];
                    sprintf(warn_buf, "Only one electric field term supported for each dimension");
                    warning_error(wi, warn_buf);
                }
                omega = std::strtod(sx[1].c_str(), &pos);
                sigma = std::strtod(sx[2].c_str(), &pos);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

            setFieldTerm(dim, aa, omega, t0, sigma);
        }
        else
        {
            setFieldTerm(dim, aa, 0, t0, 0);
        }
    }
}

void ElectricField::broadCast(const t_commrec *cr)
{
    rvec a1, omega1, sigma1, t01;
    for (int m = 0; m < DIM; m++)
    {
        a1[m]     = a(m);
        omega1[m] = omega(m);
        sigma1[m] = sigma(m);
        t01[m]    = t0(m);
    }
    gmx_bcast(DIM*sizeof(a1[0]), a1, cr);
    gmx_bcast(DIM*sizeof(omega1[0]), omega1, cr);
    gmx_bcast(DIM*sizeof(t01[0]), t01, cr);
    gmx_bcast(DIM*sizeof(sigma1[0]), sigma1, cr);

    if (!MASTER(cr))
    {
        for (int m = 0; m < DIM; m++)
        {
            setFieldTerm(m, a1[m], omega1[m], t01[m], sigma1[m]);
        }
    }
}

void ElectricField::printParameters(FILE *fp, int indent)
{
    for (int m = 0; m < DIM; m++)
    {
        indent = pr_title(fp, indent, "ElectricField");
        pr_indent(fp, indent);
        fprintf(fp, "dim = %d a = %e omega = %g t0 = %g sigma = %g\n",
                m, a(m), omega(m), t0(m), sigma(m));
    }
}

void ElectricField::setFieldTerm(int dim, real a, real omega, real t0, real sigma)
{
    range_check(dim, 0, DIM);
    efield_[dim].setField(a, omega, t0, sigma);
    isSet_ = true;
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
    return (isSet_ &&
            (efield_[XX].a() != 0 ||
             efield_[YY].a() != 0 ||
             efield_[ZZ].a() != 0));
}

void ElectricField::compare(FILE                *fp,
                            const ElectricField *f2,
                            real                 ftol,
                            real                 abstol)
{
    for (int m = 0; (m < DIM); m++)
    {
        char buf[256];

        sprintf(buf, "inputrec->field[%d]", m);
        cmp_real(fp, buf, -1, a(m), f2->a(m), ftol, abstol);
        cmp_real(fp, buf, -1, omega(m), f2->omega(m), ftol, abstol);
        cmp_real(fp, buf, -1, t0(m), f2->t0(m), ftol, abstol);
        cmp_real(fp, buf, -1, sigma(m), f2->sigma(m), ftol, abstol);
    }
}

void ElectricField::printComponents(FILE *fp, double t) const
{
    try
    {
        fprintf(fp, "%10g  %10g  %10g  %10g #FIELD\n", t,
                field(XX, t), field(YY, t), field(ZZ, t));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void calc_f_el(FILE *fplog,
               int  start, int homenr,
               real charge[], rvec f[],
               const ElectricField *efield,
               double t)
{
    if (efield->applyField())
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
        if (nullptr != fplog)
        {
            efield->printComponents(fplog, t);
        }
    }
}
