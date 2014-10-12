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
/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_H
#define GMX_EWALD_PME_H

#include <stdio.h>

#include "gromacs/legacyheaders/types/commrec_fwd.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/nrnb.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

enum {
    GMX_SUM_GRID_FORWARD, GMX_SUM_GRID_BACKWARD
};

/*! \brief Initialize \p pmedata
 *
 * Return value 0 indicates all well, non zero is an error code.
 */
int gmx_pme_init(struct gmx_pme_t **pmedata, t_commrec *cr,
                 int nnodes_major, int nnodes_minor,
                 t_inputrec *ir, int homenr,
                 gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                 gmx_bool bReproducible, int nthread);

/*! \brief Destroy the pme data structures resepectively.
 *
 * \return 0 indicates all well, non zero is an error code.
 */
int gmx_pme_destroy(FILE *log, struct gmx_pme_t **pmedata);

//@{
/*! \brief Flag values that control what gmx_pme_do() will calculate
 *
 * These can be combined with bitwise-OR if more than one thing is required.
 */
#define GMX_PME_SPREAD        (1<<0)
#define GMX_PME_SOLVE         (1<<1)
#define GMX_PME_CALC_F        (1<<2)
#define GMX_PME_CALC_ENER_VIR (1<<3)
/* This forces the grid to be backtransformed even without GMX_PME_CALC_F */
#define GMX_PME_CALC_POT      (1<<4)

/* These values label bits used for sending messages to PME nodes using the
 * routines in pme_pp.c and shouldn't conflict with the flags used there
 */
#define GMX_PME_DO_COULOMB    (1<<13)
#define GMX_PME_DO_LJ         (1<<14)

#define GMX_PME_DO_ALL_F  (GMX_PME_SPREAD | GMX_PME_SOLVE | GMX_PME_CALC_F)
//@}

/*! \brief Do a PME calculation for the long range electrostatics and/or LJ.
 *
 * The meaning of \p flags is defined above, and determines which
 * parts of the calculation are performed.
 *
 * \return 0 indicates all well, non zero is an error code.
 */
int gmx_pme_do(struct gmx_pme_t *pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real chargeA[],  real chargeB[],
               real c6A[],      real c6B[],
               real sigmaA[],   real sigmaB[],
               matrix box,      t_commrec *cr,
               int  maxshift_x, int maxshift_y,
               t_nrnb *nrnb,    gmx_wallcycle_t wcycle,
               matrix vir_q,    real ewaldcoeff_q,
               matrix vir_lj,   real ewaldcoeff_lj,
               real *energy_q,  real *energy_lj,
               real lambda_q,   real lambda_lj,
               real *dvdlambda_q, real *dvdlambda_lj,
               int flags);

/*! \brief Called on the nodes that do PME exclusively (as slaves) */
int gmx_pmeonly(struct gmx_pme_t *pme,
                t_commrec *cr,     t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                real ewaldcoeff_q, real ewaldcoeff_lj,
                t_inputrec *ir);

/*! \brief Calculate the PME grid energy V for n charges.
 *
 * The potential (found in \p pme) must have been found already with a
 * call to gmx_pme_do() with at least GMX_PME_SPREAD and GMX_PME_SOLVE
 * specified. Note that the charges are not spread on the grid in the
 * pme struct. Currently does not work in parallel or with free
 * energy.
 */
void gmx_pme_calc_energy(struct gmx_pme_t *pme, int n, rvec *x, real *q, real *V);

/*! \brief Send the charges and maxshift to out PME-only node. */
void gmx_pme_send_parameters(t_commrec *cr,
                             const interaction_const_t *ic,
                             gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                             real *chargeA, real *chargeB,
                             real *sqrt_c6A, real *sqrt_c6B,
                             real *sigmaA, real *sigmaB,
                             int maxshift_x, int maxshift_y);

/*! \brief Send the coordinates to our PME-only node and request a PME calculation */
void gmx_pme_send_coordinates(t_commrec *cr, matrix box, rvec *x,
                              gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                              real lambda_q, real lambda_lj,
                              gmx_bool bEnerVir, int pme_flags,
                              gmx_int64_t step);

/*! \brief Tell our PME-only node to finish */
void gmx_pme_send_finish(t_commrec *cr);

/*! \brief Tell our PME-only node to reset all cycle and flop counters */
void gmx_pme_send_resetcounters(t_commrec *cr, gmx_int64_t step);

/*! \brief PP nodes receive the long range forces from the PME nodes */
void gmx_pme_receive_f(t_commrec *cr,
                       rvec f[], matrix vir_q, real *energy_q,
                       matrix vir_lj, real *energy_lj,
                       real *dvdlambda_q, real *dvdlambda_lj,
                       float *pme_cycles);

#endif
