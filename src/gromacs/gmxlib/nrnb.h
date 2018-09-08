/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018, by the GROMACS development team, led by
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

#ifndef GMX_GMXLIB_NRNB_H
#define GMX_GMXLIB_NRNB_H

#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"

#define eNR_NBKERNEL_NONE (-1)

enum
{
    eNR_NBKERNEL_VDW_VF,
    eNR_NBKERNEL_VDW_F,
    eNR_NBKERNEL_ELEC_VF,
    eNR_NBKERNEL_ELEC_F,
    eNR_NBKERNEL_ELEC_W3_VF,
    eNR_NBKERNEL_ELEC_W3_F,
    eNR_NBKERNEL_ELEC_W3W3_VF,
    eNR_NBKERNEL_ELEC_W3W3_F,
    eNR_NBKERNEL_ELEC_W4_VF,
    eNR_NBKERNEL_ELEC_W4_F,
    eNR_NBKERNEL_ELEC_W4W4_VF,
    eNR_NBKERNEL_ELEC_W4W4_F,
    eNR_NBKERNEL_ELEC_VDW_VF,
    eNR_NBKERNEL_ELEC_VDW_F,
    eNR_NBKERNEL_ELEC_VDW_W3_VF,
    eNR_NBKERNEL_ELEC_VDW_W3_F,
    eNR_NBKERNEL_ELEC_VDW_W3W3_VF,
    eNR_NBKERNEL_ELEC_VDW_W3W3_F,
    eNR_NBKERNEL_ELEC_VDW_W4_VF,
    eNR_NBKERNEL_ELEC_VDW_W4_F,
    eNR_NBKERNEL_ELEC_VDW_W4W4_VF,
    eNR_NBKERNEL_ELEC_VDW_W4W4_F,

    eNR_NBKERNEL_NR,                        /* Total number of interaction-specific kernel entries */

    eNR_NBKERNEL_GENERIC = eNR_NBKERNEL_NR, /* Reuse number; KERNEL_NR is not an entry itself */
    eNR_NBKERNEL_GENERIC_CG,
    eNR_NBKERNEL_FREE_ENERGY,               /* Add other generic kernels _before_ the free energy one */

    eNR_NBKERNEL_TOTAL_NR,
    eNR_NBKERNEL_ALLVSALL = eNR_NBKERNEL_TOTAL_NR, // Reuse the symbolic constant that indicates the last kernel

    eNR_NBNXN_DIST2,
    eNR_NBNXN_LJ_RF,          eNR_NBNXN_LJ_RF_E,
    eNR_NBNXN_LJ_TAB,         eNR_NBNXN_LJ_TAB_E,
    eNR_NBNXN_LJ_EWALD,       eNR_NBNXN_LJ_EWALD_E,
    eNR_NBNXN_LJ,             eNR_NBNXN_LJ_E,
    eNR_NBNXN_RF,             eNR_NBNXN_RF_E,
    eNR_NBNXN_TAB,            eNR_NBNXN_TAB_E,
    eNR_NBNXN_EWALD,          eNR_NBNXN_EWALD_E,
    eNR_NBNXN_ADD_LJ_FSW,     eNR_NBNXN_ADD_LJ_FSW_E,
    eNR_NBNXN_ADD_LJ_PSW,     eNR_NBNXN_ADD_LJ_PSW_E,
    eNR_NBNXN_ADD_LJ_EWALD,   eNR_NBNXN_ADD_LJ_EWALD_E,
    eNR_NB14,
    eNR_WEIGHTS,              eNR_SPREAD,               eNR_SPREADBSP,
    eNR_GATHERF,              eNR_GATHERFBSP,           eNR_FFT,
    eNR_CONV,                 eNR_SOLVEPME, eNR_NS,      eNR_RESETX,
    eNR_SHIFTX,               eNR_CGCM,                 eNR_FSUM,
    eNR_BONDS,                eNR_G96BONDS,             eNR_FENEBONDS,
    eNR_TABBONDS,             eNR_RESTRBONDS,           eNR_LINEAR_ANGLES,
    eNR_ANGLES,               eNR_G96ANGLES,            eNR_QANGLES,
    eNR_TABANGLES,            eNR_PROPER,               eNR_IMPROPER,
    eNR_RB,                   eNR_FOURDIH,              eNR_TABDIHS,
    eNR_DISRES,               eNR_ORIRES,               eNR_DIHRES,
    eNR_POSRES,               eNR_FBPOSRES,
    eNR_ANGRES,               eNR_ANGRESZ,
    eNR_MORSE,                eNR_CUBICBONDS,           eNR_WALLS,
    eNR_POLARIZE,             eNR_ANHARM_POL,
    eNR_WPOL,                 eNR_THOLE,                eNR_VIRIAL,
    eNR_UPDATE,               eNR_EXTUPDATE,            eNR_STOPCM,
    eNR_PCOUPL,               eNR_EKIN,                 eNR_LINCS,
    eNR_LINCSMAT,             eNR_SHAKE,                eNR_CONSTR_V,
    eNR_SHAKE_RIJ,            eNR_CONSTR_VIR,           eNR_SETTLE,
    eNR_VSITE2,               eNR_VSITE3,               eNR_VSITE3FD,
    eNR_VSITE3FAD,            eNR_VSITE3OUT,            eNR_VSITE4FD,
    eNR_VSITE4FDN,            eNR_VSITEN,
    eNR_CMAP,                 eNR_UREY_BRADLEY,         eNR_CROSS_BOND_BOND,
    eNR_CROSS_BOND_ANGLE,
    eNRNB
};


typedef struct t_nrnb
{
    double n[eNRNB];
}
t_nrnb;
struct t_commrec;
struct t_nrnb;

void init_nrnb(t_nrnb *nrnb);

void cp_nrnb(t_nrnb *dest, t_nrnb *src);

void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2);

void print_nrnb(FILE *out, t_nrnb *nrnb);

void _inc_nrnb(t_nrnb *nrnb, int enr, int inc, char *file, int line);

#ifdef DEBUG_NRNB
#define inc_nrnb(nrnb, enr, inc) _inc_nrnb(nrnb, enr, inc, __FILE__, __LINE__)
#else
#define inc_nrnb(nrnb, enr, inc) (nrnb)->n[enr] += inc
#endif


void print_flop(FILE *out, t_nrnb *nrnb, double *nbfs, double *mflop);
/* Calculates the non-bonded forces and flop count.
 * When out!=NULL also prints the full count table.
 */

void print_perf(FILE *out, double nodetime, double realtime,
                int64_t nsteps, double delta_t,
                double nbfs, double mflop);
/* Prints the performance, nbfs and mflop come from print_flop */

void pr_load(FILE *log, struct t_commrec *cr, t_nrnb nrnb[]);
/* Print detailed load balancing info */

int cost_nrnb(int enr);
/* Cost in i860 cycles of this component of MD */

const char *nrnb_str(int enr);
/* Name of this component */

#endif  /* GMX_GMXLIB_NRNB_H */
