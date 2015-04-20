/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef GMX_LEGACYHEADERS_TYPES_QMMMREC_H
#define GMX_LEGACYHEADERS_TYPES_QMMMREC_H

#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int                nrQMatoms;      /* total nr of QM atoms              */
    rvec              *xQM;            /* shifted to center of box          */
    int               *indexQM;        /* atom i = atom indexQM[i] in mdrun */
    int               *atomicnumberQM; /* atomic numbers of QM atoms        */
    real              *QMcharges;      /* atomic charges of QM atoms(ONIOM) */
    int               *shiftQM;
    int                QMcharge;       /* charge of the QM system           */
    int                multiplicity;   /* multipicity (no of unpaired eln)  */
    int                QMmethod;       /* see enums.h for all methods       */
    int                QMbasis;        /* see enums.h for all bases         */
    int                nelectrons;     /* total number of elecs in QM region*/
    gmx_bool           bTS;            /* Optimize a TS, only steep, no md  */
    gmx_bool           bOPT;           /* Optimize QM subsys, only steep, no md  */
    gmx_bool          *frontatoms;     /* qm atoms on the QM side of a QM-MM bond */
    /* Gaussian specific stuff */
    int                nQMcpus;        /* no. of CPUs used for the QM calc. */
    int                QMmem;          /* memory for the gaussian calc.     */
    int                accuracy;       /* convergence criterium (E(-x))     */
    gmx_bool           cpmcscf;        /* using cpmcscf(l1003)*/
    char              *gauss_dir;
    char              *gauss_exe;
    char              *devel_dir;
    char              *orca_basename; /* basename for I/O with orca        */
    char              *orca_dir;      /* directory for ORCA                */
    real              *c6;
    real              *c12;
    /* Surface hopping stuff */
    gmx_bool           bSH;     /* surface hopping (diabatic only)   */
    real               SAon;    /* at which energy gap the SA starts */
    real               SAoff;   /* at which energy gap the SA stops  */
    int                SAsteps; /* stepwise switchinng on the SA     */
    int                SAstep;  /* current state of SA               */
    int                CIdim;
    real              *CIvec1;
    real              *CIvec2;
    real              *CIvec1old;
    real              *CIvec2old;
    ivec               SHbasis;
    int                CASelectrons;
    int                CASorbitals;
} t_QMrec;

typedef struct {
    int            nrMMatoms;   /* nr of MM atoms, updated every step*/
    rvec          *xMM;         /* shifted to center of box          */
    int           *indexMM;     /* atom i = atom indexMM[I] in mdrun */
    real          *MMcharges;   /* MM point charges in std QMMM calc.*/
    int           *shiftMM;
    int           *MMatomtype;  /* only important for semi-emp.      */
    real           scalefactor;
    /* gaussian specific stuff */
    real          *c6;
    real          *c12;
} t_MMrec;


typedef struct {
    int             QMMMscheme; /* ONIOM (multi-layer) or normal          */
    int             nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
    t_QMrec       **qm;         /* atoms and run params for each QM group */
    t_MMrec        *mm;         /* there can only be one MM subsystem !   */
} t_QMMMrec;

#ifdef __cplusplus
}
#endif
#endif
