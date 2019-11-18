/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "qm_gamess.h"

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

/* QMMM sub routines */
/* mopac interface routines */


static void F77_FUNC(inigms, IMIGMS)();

static void F77_FUNC(grads, GRADS)(const int*  nrqmat,
                                   real*       qmcrd,
                                   const int*  nrmmat,
                                   const real* mmchrg,
                                   real*       mmcrd,
                                   real*       qmgrad,
                                   real*       mmgrad,
                                   real*       energy);

#if !GMX_QMMM_GAMESS
// Stub definitions to make compilation succeed when not configured
// for GAMESS support. In that case, the module gives a fatal error
// when the initialization function is called, so there is no need to
// issue fatal errors here, because that introduces problems with
// tools suggesting and prohibiting noreturn attributes.

void F77_FUNC(inigms, IMIGMS)(){};
// NOLINTNEXTLINE(readability-named-parameter)
void F77_FUNC(grads, GRADS)(const int*, real*, const int*, const real*, real*, real*, real*, real*){};
#endif

void init_gamess(const t_commrec* cr, t_QMrec* qm, t_MMrec* mm)
{
    /* it works hopelessly complicated :-)
     * first a file is written. Then the standard gamess input/output
     * routine is called (no system()!) to set up all fortran arrays.
     * this routine writes a punch file, like in a normal gamess run.
     * via this punch file the other games routines, needed for gradient
     * and energy evaluations are called. This setup works fine for
     * dynamics simulations. 7-6-2002 (London)
     */
    int   i, j;
    FILE* out;
    char  periodic_system[37][3] = { "XX", "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ",
                                    "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ",
                                    "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                                    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr" };
    if (!GMX_QMMM_GAMESS)
    {
        gmx_fatal(FARGS,
                  "Cannot call GAMESS unless linked against it. Use cmake "
                  "-DGMX_QMMM_PROGRAM=GAMESS, and ensure that linking will work correctly.");
    }

    if (PAR(cr))
    {

        if (MASTER(cr))
        {
            out = fopen("FOR009", "w");
            /* of these options I am not completely sure....  the overall
             * preformance on more than 4 cpu's is rather poor at the moment.
             */
            fprintf(out, "memory 48000000\nPARALLEL IOMODE SCREENED\n");
            fprintf(out, "ELEC %d\nMULT %d\nSUPER ON\nNOSYM\nGEOMETRY ANGSTROM\n", qm->nelectrons,
                    qm->multiplicity);
            for (i = 0; i < qm->nrQMatoms; i++)
            {
#ifdef DOUBLE
                fprintf(out, "%10.7lf  %10.7lf  %10.7lf  %5.3lf  %2s\n", i / 2., i / 3., i / 4.,
                        qm->atomicnumberQM[i] * 1.0, periodic_system[qm->atomicnumberQM[i]]);
#else
                fprintf(out, "%10.7f  %10.7f  %10.7f  %5.3f  %2s\n", i / 2., i / 3., i / 4.,
                        qm->atomicnumberQM[i] * 1.0, periodic_system[qm->atomicnumberQM[i]]);
#endif
            }
            if (mm->nrMMatoms)
            {
                for (j = i; j < i + 2; j++)
                {
#ifdef DOUBLE
                    fprintf(out, "%10.7lf  %10.7lf  %10.7lf  %5.3lf  BQ\n", j / 5., j / 6., j / 7., 1.0);
#else
                    fprintf(out, "%10.7f  %10.7f  %10.7f  %5.3f  BQ\n", j / 5., j / 6., j / 7., 2.0);
#endif
                }
            }
            fprintf(out, "END\nBASIS %s\nRUNTYPE GRADIENT\nSCFTYPE %s\n",
                    eQMbasis_names[qm->QMbasis], eQMmethod_names[qm->QMmethod]); /* see enum.h */
            fclose(out);
        }
        gmx_barrier(cr);
        F77_FUNC(inigms, IMIGMS)();
    }
    else /* normal serial run */

    {
        out = fopen("FOR009", "w");
        /* of these options I am not completely sure....  the overall
         * preformance on more than 4 cpu's is rather poor at the moment.
         */
        fprintf(out, "ELEC %d\nMULT %d\nSUPER ON\nNOSYM\nGEOMETRY ANGSTROM\n", qm->nelectrons,
                qm->multiplicity);
        for (i = 0; i < qm->nrQMatoms; i++)
        {
#ifdef DOUBLE
            fprintf(out, "%10.7lf  %10.7lf  %10.7lf  %5.3lf  %2s\n", i / 2., i / 3., i / 4.,
                    qm->atomicnumberQM[i] * 1.0, periodic_system[qm->atomicnumberQM[i]]);
#else
            fprintf(out, "%10.7f  %10.7f  %10.7f  %5.3f  %2s\n", i / 2., i / 3., i / 4.,
                    qm->atomicnumberQM[i] * 1.0, periodic_system[qm->atomicnumberQM[i]]);
#endif
        }
        if (mm->nrMMatoms)
        {
            for (j = i; j < i + 2; j++)
            {
#ifdef DOUBLE
                fprintf(out, "%10.7lf  %10.7lf  %10.7lf  %5.3lf  BQ\n", j / 5., j / 6., j / 7., 1.0);
#else
                fprintf(out, "%10.7f  %10.7f  %10.7f  %5.3f  BQ\n", j / 5., j / 6., j / 7., 2.0);
#endif
            }
        }
        fprintf(out, "END\nBASIS %s\nRUNTYPE GRADIENT\nSCFTYPE %s\n", eQMbasis_names[qm->QMbasis],
                eQMmethod_names[qm->QMmethod]); /* see enum.h */
        F77_FUNC(inigms, IMIGMS)();
    }
}

real call_gamess(const t_QMrec* qm, const t_MMrec* mm, rvec f[], rvec fshift[])
{
    /* do the actual QMMM calculation using GAMESS-UK. In this
     * implementation (3-2001) a system call is made to the GAMESS-UK
     * binary. Now we are working to get the electron integral, SCF, and
     * gradient routines linked directly
     */
    int  i, j;
    real QMener = 0.0, *qmgrad, *mmgrad, *mmcrd, *qmcrd, energy = 0;

    snew(qmcrd, 3 * (qm->nrQMatoms));
    snew(mmcrd, 3 * (mm->nrMMatoms));
    snew(qmgrad, 3 * (qm->nrQMatoms));
    snew(mmgrad, 3 * (mm->nrMMatoms));

    /* copy the data from qr into the arrays that are going to be used
     * in the fortran routines of gamess
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            qmcrd[DIM * i + j] = 1 / BOHR2NM * qm->xQM[i][j];
        }
    }
    for (i = 0; i < mm->nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            mmcrd[DIM * i + j] = 1 / BOHR2NM * mm->xMM[i][j];
        }
    }
    for (i = 0; i < 3 * qm->nrQMatoms; i += 3)
    {
        fprintf(stderr, "%8.5f, %8.5f, %8.5f\n", qmcrd[i], qmcrd[i + 1], qmcrd[i + 2]);
    }

    F77_FUNC(grads, GRADS)
    (&qm->nrQMatoms, qmcrd, &mm->nrMMatoms, mm->MMcharges, mmcrd, qmgrad, mmgrad, &energy);

    for (i = 0; i < qm->nrQMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD * qmgrad[3 * i + j];
            fshift[i][j] = HARTREE_BOHR2MD * qmgrad[3 * i + j];
        }
    }
    for (i = 0; i < mm->nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD * mmgrad[3 * i + j];
            fshift[i][j] = HARTREE_BOHR2MD * mmgrad[3 * i + j];
        }
    }
    /* convert a.u to kJ/mol */
    QMener = energy * HARTREE2KJ * AVOGADRO;
    return (QMener);
}

#pragma GCC diagnostic pop
