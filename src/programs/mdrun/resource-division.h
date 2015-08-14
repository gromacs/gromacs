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

#ifndef GMX_RESOURCE_DIVISION_H
#define GMX_RESOURCE_DIVISION_H

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec_fwd.h"

/* Return the number of threads to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than #GPU.
 * If necessary, this function will modify hw_opt->nthreads_omp.
 */
int get_nthreads_mpi(const gmx_hw_info_t *hwinfo,
                     gmx_hw_opt_t        *hw_opt,
                     const t_inputrec    *inputrec,
                     const gmx_mtop_t    *mtop,
                     const t_commrec     *cr,
                     FILE                *fplog);

/* Check if the number of OpenMP threads is within reasonable range
 * considering the hardware used. This is a crude check, but mainly
 * intended to catch cases where the user starts 1 MPI rank per hardware
 * thread or 1 rank per physical node.
 * With a sub-optimal setup a note is printed to fplog and stderr when
 * bNtOmpSet==TRUE; with bNtOptOptionSet==FALSE a fatal error is issued.
 * This function should be called after thread-MPI and OpenMP are set up.
 */
void check_resource_division_efficiency(const gmx_hw_info_t *hwinfo,
                                        const gmx_hw_opt_t  *hw_opt,
                                        gmx_bool             bNtOmpOptionSet,
                                        t_commrec           *cr,
                                        FILE                *fplog);

/* Checks we can do when we don't (yet) know the cut-off scheme */
void check_and_update_hw_opt_1(gmx_hw_opt_t    *hw_opt,
                               const t_commrec *cr);

/* Checks we can do when we know the cut-off scheme */
void check_and_update_hw_opt_2(gmx_hw_opt_t *hw_opt,
                               int           cutoff_scheme);

/* Checks we can do when we know the thread-MPI rank count */
void check_and_update_hw_opt_3(gmx_hw_opt_t *hw_opt);

#endif /* GMX_RESOURCE_DIVISION_H */
