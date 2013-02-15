/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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


#ifdef __cplusplus
extern "C" {
#endif

/* define USE_MPE if you want MPE logging
 *
 * you then need to link with the appropriate libraries
 * that come with the mpich distribution (can be found in the
 * mpe subdirectory */
/* #define USE_MPE */

/* define BARRIERS if you want to have extra MPI_Barriers
 * in the code which might help analyzing the MPE logfiles
 */
/* #define BARRIERS */
#ifdef BARRIERS
#define GMX_BARRIER(communicator) MPI_Barrier(communicator)
#else
#define GMX_BARRIER(communicator)
#endif

#ifdef USE_MPE
#define GMX_MPE_LOG(event) MPE_Log_event(event, 0, "")
#else
#define GMX_MPE_LOG(event)
#endif

#ifdef USE_MPE
#include <mpe.h>
/* Define MPE logging events here */
/* General events */
int ev_timestep1,               ev_timestep2;
int ev_ns_start,                ev_ns_finish;
int ev_calc_bonds_start,        ev_calc_bonds_finish;
int ev_send_coordinates_start,  ev_send_coordinates_finish;
int ev_update_fr_start,         ev_update_fr_finish;
int ev_clear_rvecs_start,       ev_clear_rvecs_finish;
int ev_output_start,            ev_output_finish;
int ev_update_start,            ev_update_finish;
int ev_force_start,             ev_force_finish;
int ev_do_fnbf_start,           ev_do_fnbf_finish;

/* Shift related events*/
int ev_shift_start,             ev_shift_finish;
int ev_unshift_start,           ev_unshift_finish;
int ev_mk_mshift_start,         ev_mk_mshift_finish;

/* PME related events */
int ev_pme_start,               ev_pme_finish;
int ev_spread_on_grid_start,    ev_spread_on_grid_finish;
int ev_sum_qgrid_start,         ev_sum_qgrid_finish;
int ev_gmxfft3d_start,          ev_gmxfft3d_finish;
int ev_solve_pme_start,         ev_solve_pme_finish;
int ev_gather_f_bsplines_start, ev_gather_f_bsplines_finish;
int ev_reduce_start,            ev_reduce_finish;
int ev_rscatter_start,          ev_rscatter_finish;
int ev_alltoall_start,          ev_alltoall_finish;
int ev_pmeredist_start,         ev_pmeredist_finish;
int ev_init_pme_start,          ev_init_pme_finish;
int ev_global_stat_start,       ev_global_stat_finish;
int ev_sum_lrforces_start,      ev_sum_lrforces_finish;
int ev_virial_start,            ev_virial_finish;
int ev_sort_start,              ev_sort_finish;
int ev_sum_qgrid_start,         ev_sum_qgrid_finish;

/* Essential dynamics related events */
int ev_edsam_start,             ev_edsam_finish;
int ev_get_group_x_start,       ev_get_group_x_finish;
int ev_ed_apply_cons_start,     ev_ed_apply_cons_finish;
int ev_fit_to_reference_start,  ev_fit_to_reference_finish;
#endif

#ifdef __cplusplus
}
#endif
