/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _tgroup_h
#define _tgroup_h

#include "typedefs.h"
#include "network.h"

#ifdef __cplusplus
extern "C" {
#endif

void
init_temperature_coupling_outputs(const t_grpopts *opts,
                                  gmx_temperature_coupling_outputs_t *temperature_coupling_outputs);
/* Allocate memory and set up. */

void
init_constant_acceleration(const gmx_mtop_t *mtop,
                           const t_grpopts *opts,
                           gmx_constant_acceleration_t *const_acc);
/* Allocate memory and set up. */

void init_ekindata(const t_grpopts *opts,
                   gmx_ekindata_t  *ekind);
/* Allocate memory and set up. */

void done_ekindata(gmx_ekindata_t *ekind);
/* Free the memory */

void accumulate_u(t_commrec *cr,
                  gmx_constant_acceleration_t *constant_acceleration);
/* Accumulate mean velocities of particles over all ranks?
 * TODO fold this into the normal energy communication (if it survives) */

real sum_ekin(t_grpopts *opts,
              gmx_ekindata_t *ekind,
              gmx_temperature_coupling_outputs_t *temperature_coupling_outputs,
              real *dekindlambda,
              gmx_bool bEkinFullStep, gmx_bool bScaleEkin);
/* Sum the group ekins into total ekin and calc temp per group,
 * return total temperature.
 */

void update_ekindata(int start, int homenr,
                     gmx_temperature_coupling_outputs_t *temperature_coupling_outputs,
                     gmx_constant_acceleration_t *constant_acceleration,
                     t_grpopts *opts, rvec v[], t_mdatoms *md, real lambda);
/* Do the update of group velocities (if constant acceleration) and
 * (partial) group ekin.
 */

#ifdef __cplusplus
}
#endif

#endif  /* _tgroup_h */
