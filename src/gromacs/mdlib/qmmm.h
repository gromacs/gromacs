/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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
#ifndef GMX_MDLIB_QMMM_H
#define GMX_MDLIB_QMMM_H

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/tgroup.h"

struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_QMMMrec;

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
    /* Surface hopping stuff */
    gmx_bool           bSH;           /* surface hopping (diabatic only)   */
    real               SAon;          /* at which energy gap the SA starts */
    real               SAoff;         /* at which energy gap the SA stops  */
    int                SAsteps;       /* stepwise switchinng on the SA     */
    int                SAstep;        /* current state of SA               */
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
} t_MMrec;


typedef struct t_QMMMrec {
    int             QMMMscheme; /* ONIOM (multi-layer) or normal          */
    int             nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
    t_QMrec       **qm;         /* atoms and run params for each QM group */
    t_MMrec        *mm;         /* there can only be one MM subsystem !   */
} t_QMMMrec;

void atomic_number(int nr, char ***atomtype, int *nucnum);

t_QMMMrec *mk_QMMMrec(void);
/* allocates memory for QMMMrec */

void init_QMMMrec(t_commrec  *cr,
                  gmx_mtop_t *mtop,
                  t_inputrec *ir,
                  t_forcerec *fr);

/* init_QMMMrec initializes the QMMM record. From
 * topology->atoms.atomname and topology->atoms.atomtype the atom
 * names and types are read; from inputrec->QMcharge
 * resp. inputrec->QMmult the nelecs and multiplicity are determined
 * and md->cQMMM gives numbers of the MM and QM atoms
 */

void update_QMMMrec(t_commrec       *cr,
                    t_forcerec      *fr,
                    const rvec      *x,
                    const t_mdatoms *md,
                    const matrix     box);

/* update_QMMMrec fills the MM stuff in QMMMrec. The MM atoms are
 * taken froom the neighbourlists of the QM atoms. In a QMMM run this
 * routine should be called at every step, since it updates the MM
 * elements of the t_QMMMrec struct.
 */

real calculate_QMMM(t_commrec  *cr,
                    rvec        f[],
                    t_forcerec *fr);

/* QMMM computes the QM forces. This routine makes either function
 * calls to gmx QM routines (derived from MOPAC7 (semi-emp.) and MPQC
 * (ab initio)) or generates input files for an external QM package
 * (listed in QMMMrec.QMpackage). The binary of the QM package is
 * called by system().
 */

#endif
