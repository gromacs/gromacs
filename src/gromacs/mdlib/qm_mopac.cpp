/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "config.h"

#if GMX_QMMM_MOPAC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"


/* mopac interface routines */
void
    F77_FUNC(domldt, DOMLDT) (int *nrqmat, int labels[], char keywords[]);

void
    F77_FUNC(domop, DOMOP) (int *nrqmat, double qmcrd[], int *nrmmat,
                            double mmchrg[], double mmcrd[], double qmgrad[],
                            double mmgrad[], double *energy, double qmcharges[]);



void init_mopac(t_QMrec *qm)
{
    /* initializes the mopac routines ans sets up the semiempirical
     * computation by calling moldat(). The inline mopac routines can
     * only perform gradient operations. If one would like to optimize a
     * structure or find a transition state at PM3 level, gaussian is
     * used instead.
     */
    char
    *keywords;

    snew(keywords, 240);

    if (!qm->bSH)  /* if rerun then grad should not be done! */
    {
        sprintf(keywords, "PRECISE GEO-OK CHARGE=%d GRAD MMOK ANALYT %s\n",
                qm->QMcharge,
                eQMmethod_names[qm->QMmethod]);
    }
    else
    {
        sprintf(keywords, "PRECISE GEO-OK CHARGE=%d SINGLET GRAD %s C.I.=(%d,%d) root=2 MECI \n",
                qm->QMcharge,
                eQMmethod_names[qm->QMmethod],
                qm->CASorbitals, qm->CASelectrons/2);
    }
    F77_FUNC(domldt, DOMLDT) (&qm->nrQMatoms, qm->atomicnumberQM, keywords);
    fprintf(stderr, "keywords are: %s\n", keywords);
    free(keywords);

} /* init_mopac */

real call_mopac(t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
    /* do the actual QMMM calculation using directly linked mopac subroutines
     */
    double /* always double as the MOPAC routines are always compiled in
              double precission! */
    *qmcrd = NULL, *qmchrg = NULL, *mmcrd = NULL, *mmchrg = NULL,
    *qmgrad, *mmgrad = NULL, energy;
    int
        i, j;
    real
        QMener = 0.0;
    snew(qmcrd, 3*(qm->nrQMatoms));
    snew(qmgrad, 3*(qm->nrQMatoms));
    /* copy the data from qr into the arrays that are going to be used
     * in the fortran routines of MOPAC
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            qmcrd[3*i+j] = (double)qm->xQM[i][j]*10;
        }
    }
    if (mm->nrMMatoms)
    {
        /* later we will add the point charges here. There are some
         * conceptual problems with semi-empirical QM in combination with
         * point charges that we need to solve first....
         */
        gmx_fatal(FARGS, "At present only ONIOM is allowed in combination"
                  " with MOPAC QM subroutines\n");
    }
    else
    {
        /* now compute the energy and the gradients.
         */

        snew(qmchrg, qm->nrQMatoms);
        F77_FUNC(domop, DOMOP) (&qm->nrQMatoms, qmcrd, &mm->nrMMatoms,
                                mmchrg, mmcrd, qmgrad, mmgrad, &energy, qmchrg);
        /* add the gradients to the f[] array, and also to the fshift[].
         * the mopac gradients are in kCal/angstrom.
         */
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                f[i][j]       = (real)10*CAL2JOULE*qmgrad[3*i+j];
                fshift[i][j]  = (real)10*CAL2JOULE*qmgrad[3*i+j];
            }
        }
        QMener = (real)CAL2JOULE*energy;
        /* do we do something with the mulliken charges?? */

        free(qmchrg);
    }
    free(qmgrad);
    free(qmcrd);
    return (QMener);
}

real call_mopac_SH(t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
    /* do the actual SH QMMM calculation using directly linked mopac
       subroutines */

    double /* always double as the MOPAC routines are always compiled in
              double precission! */
    *qmcrd = NULL, *qmchrg = NULL, *mmcrd = NULL, *mmchrg = NULL,
    *qmgrad, *mmgrad = NULL, energy;
    int
        i, j;
    real
        QMener = 0.0;

    snew(qmcrd, 3*(qm->nrQMatoms));
    snew(qmgrad, 3*(qm->nrQMatoms));
    /* copy the data from qr into the arrays that are going to be used
     * in the fortran routines of MOPAC
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            qmcrd[3*i+j] = (double)qm->xQM[i][j]*10;
        }
    }
    if (mm->nrMMatoms)
    {
        /* later we will add the point charges here. There are some
         * conceptual problems with semi-empirical QM in combination with
         * point charges that we need to solve first....
         */
        gmx_fatal(FARGS, "At present only ONIOM is allowed in combination with MOPAC\n");
    }
    else
    {
        /* now compute the energy and the gradients.
         */
        snew(qmchrg, qm->nrQMatoms);

        F77_FUNC(domop, DOMOP) (&qm->nrQMatoms, qmcrd, &mm->nrMMatoms,
                                mmchrg, mmcrd, qmgrad, mmgrad, &energy, qmchrg);
        /* add the gradients to the f[] array, and also to the fshift[].
         * the mopac gradients are in kCal/angstrom.
         */
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                f[i][j]      = (real)10*CAL2JOULE*qmgrad[3*i+j];
                fshift[i][j] = (real)10*CAL2JOULE*qmgrad[3*i+j];
            }
        }
        QMener = (real)CAL2JOULE*energy;
    }
    free(qmgrad);
    free(qmcrd);
    return (QMener);
} /* call_mopac_SH */

#else
int
    gmx_qmmm_mopac_empty;
#endif
