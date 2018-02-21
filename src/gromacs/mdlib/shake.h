/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#ifndef GMX_MDLIB_SHAKE_H
#define GMX_MDLIB_SHAKE_H

#include "gromacs/mdlib/constr.h"
#include "gromacs/topology/idef.h"

struct t_inputrec;

/* Abstract type for SHAKE that is defined only in the file that uses it */
typedef struct gmx_shakedata *gmx_shakedata_t;

gmx_shakedata_t shake_init();
/* Initializes and return the SHAKE data structure */

gmx_bool bshakef(FILE           *log,          /* Log file			*/
                 gmx_shakedata_t shaked,       /* Total number of atoms	*/
                 real            invmass[],    /* Atomic masses		*/
                 int             nblocks,      /* The number of shake blocks	*/
                 int             sblock[],     /* The shake blocks             */
                 t_idef         *idef,         /* The interaction def		*/
                 t_inputrec     *ir,           /* Input record		        */
                 rvec            x_s[],        /* Coords before update		*/
                 rvec            prime[],      /* Output coords		*/
                 t_nrnb         *nrnb,         /* Performance measure          */
                 real           *lagr,         /* The Lagrange multipliers     */
                 real            lambda,       /* FEP lambda                   */
                 real           *dvdlambda,    /* FEP force                    */
                 real            invdt,        /* 1/delta_t                    */
                 rvec           *v,            /* Also constrain v if v!=NULL  */
                 gmx_bool        bCalcVir,     /* Calculate r x m delta_r      */
                 tensor          vir_r_m_dr,   /* sum r x m delta_r            */
                 gmx_bool        bDumpOnError, /* Dump debugging stuff on error*/
                 int             econq);       /* which type of constraint is occurring */
/* Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return TRUE when OK, FALSE when shake-error
 */

void cshake(const int iatom[], int ncon, int *nnit, int maxnit,
            const real dist2[], real xp[], const real rij[], const real m2[], real omega,
            const real invmass[], const real tt[], real lagr[], int *nerror);
/* Regular iterative shake */

#endif
