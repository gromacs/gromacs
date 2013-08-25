/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _pme_h
#define _pme_h

#include <stdio.h>
#include "visibility.h"
#include "typedefs.h"
#include "gmxcomplex.h"
#include "gmx_wallcycle.h"
#include "sim_util.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef real *splinevec[DIM];

enum {
    GMX_SUM_QGRID_FORWARD, GMX_SUM_QGRID_BACKWARD
};

GMX_LIBMD_EXPORT
int gmx_pme_init(gmx_pme_t *pmedata, t_commrec *cr,
                 int nnodes_major, int nnodes_minor,
                 t_inputrec *ir, int homenr,
                 gmx_bool bFreeEnergy, gmx_bool bReproducible, int nthread);
/* Initialize the pme data structures resepectively.
 * Return value 0 indicates all well, non zero is an error code.
 */

GMX_LIBMD_EXPORT
int gmx_pme_reinit(gmx_pme_t *         pmedata,
                   t_commrec *         cr,
                   gmx_pme_t           pme_src,
                   const t_inputrec *  ir,
                   ivec                grid_size);
/* As gmx_pme_init, but takes most settings, except the grid, from pme_src */

int gmx_pme_destroy(FILE *log, gmx_pme_t *pmedata);
/* Destroy the pme data structures resepectively.
 * Return value 0 indicates all well, non zero is an error code.
 */

#define GMX_PME_SPREAD_Q      (1<<0)
#define GMX_PME_SOLVE         (1<<1)
#define GMX_PME_CALC_F        (1<<2)
#define GMX_PME_CALC_ENER_VIR (1<<3)
/* This forces the grid to be backtransformed even without GMX_PME_CALC_F */
#define GMX_PME_CALC_POT      (1<<4)
#define GMX_PME_DO_ALL_F  (GMX_PME_SPREAD_Q | GMX_PME_SOLVE | GMX_PME_CALC_F)

int gmx_pme_do(gmx_pme_t pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real chargeA[],  real chargeB[],
               matrix box,      t_commrec *cr,
               int  maxshift_x, int maxshift_y,
               t_nrnb *nrnb,    gmx_wallcycle_t wcycle,
               matrix lrvir,    real ewaldcoeff,
               real *energy,    real lambda,
               real *dvdlambda, int flags);
/* Do a PME calculation for the long range electrostatics.
 * flags, defined above, determine which parts of the calculation are performed.
 * Return value 0 indicates all well, non zero is an error code.
 */

GMX_LIBMD_EXPORT
int gmx_pmeonly(gmx_pme_t pme,
                t_commrec *cr,     t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                gmx_runtime_t *runtime,
                real ewaldcoeff,   gmx_bool bGatherOnly,
                t_inputrec *ir);
/* Called on the nodes that do PME exclusively (as slaves)
 */

void gmx_pme_calc_energy(gmx_pme_t pme, int n, rvec *x, real *q, real *V);
/* Calculate the PME grid energy V for n charges with a potential
 * in the pme struct determined before with a call to gmx_pme_do
 * with at least GMX_PME_SPREAD_Q and GMX_PME_SOLVE specified.
 * Note that the charges are not spread on the grid in the pme struct.
 * Currently does not work in parallel or with free energy.
 */

/* The following three routines are for PME/PP node splitting in pme_pp.c */

/* Abstract type for PME <-> PP communication */
typedef struct gmx_pme_pp *gmx_pme_pp_t;

GMX_LIBMD_EXPORT
void gmx_pme_check_restrictions(int pme_order,
                                int nkx, int nky, int nkz,
                                int nnodes_major,
                                int nnodes_minor,
                                gmx_bool bUseThreads,
                                gmx_bool bFatal,
                                gmx_bool *bValidSettings);
/* Check restrictions on pme_order and the PME grid nkx,nky,nkz.
 * With bFatal=TRUE, a fatal error is generated on violation,
 * bValidSettings=NULL can be passed.
 * With bFatal=FALSE, *bValidSettings reports the validity of the settings.
 * bUseThreads tells if any MPI rank doing PME uses more than 1 threads.
 * If at calling you bUseThreads is unknown, pass TRUE for conservative
 * checking.
 */

gmx_pme_pp_t gmx_pme_pp_init(t_commrec *cr);
/* Initialize the PME-only side of the PME <-> PP communication */

void gmx_pme_send_q(t_commrec *cr,
                    gmx_bool bFreeEnergy, real *chargeA, real *chargeB,
                    int maxshift_x, int maxshift_y);
/* Send the charges and maxshift to out PME-only node. */

void gmx_pme_send_x(t_commrec *cr, matrix box, rvec *x,
                    gmx_bool bFreeEnergy, real lambda,
                    gmx_bool bEnerVir,
                    gmx_large_int_t step);
/* Send the coordinates to our PME-only node and request a PME calculation */

GMX_LIBMD_EXPORT
void gmx_pme_send_finish(t_commrec *cr);
/* Tell our PME-only node to finish */

GMX_LIBMD_EXPORT
void gmx_pme_send_switchgrid(t_commrec *cr, ivec grid_size, real ewaldcoeff);
/* Tell our PME-only node to switch to a new grid size */

GMX_LIBMD_EXPORT
void gmx_pme_send_resetcounters(t_commrec *cr, gmx_large_int_t step);
/* Tell our PME-only node to reset all cycle and flop counters */

void gmx_pme_receive_f(t_commrec *cr,
                       rvec f[], matrix vir,
                       real *energy, real *dvdlambda,
                       float *pme_cycles);
/* PP nodes receive the long range forces from the PME nodes */

/* Return values for gmx_pme_recv_q_x */
enum {
    pmerecvqxX,            /* calculate PME mesh interactions for new x    */
    pmerecvqxFINISH,       /* the simulation should finish, we should quit */
    pmerecvqxSWITCHGRID,   /* change the PME grid size                     */
    pmerecvqxRESETCOUNTERS /* reset the cycle and flop counters            */
};

int gmx_pme_recv_q_x(gmx_pme_pp_t pme_pp,
                     int *natoms,
                     real **chargeA, real **chargeB,
                     matrix box, rvec **x, rvec **f,
                     int *maxshift_x, int *maxshift_y,
                     gmx_bool *bFreeEnergy, real *lambda,
                     gmx_bool *bEnerVir,
                     gmx_large_int_t *step,
                     ivec grid_size, real *ewaldcoeff);
;
/* With return value:
 * pmerecvqxX:             all parameters set, chargeA and chargeB can be NULL
 * pmerecvqxFINISH:        no parameters set
 * pmerecvqxSWITCHGRID:    only grid_size and *ewaldcoeff are set
 * pmerecvqxRESETCOUNTERS: *step is set
 */

void gmx_pme_send_force_vir_ener(gmx_pme_pp_t pme_pp,
                                 rvec *f, matrix vir,
                                 real energy, real dvdlambda,
                                 float cycles);
/* Send the PME mesh force, virial and energy to the PP-only nodes */

#ifdef __cplusplus
}
#endif

#endif
