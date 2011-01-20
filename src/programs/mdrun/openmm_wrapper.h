/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _OPENMM_WRAPPER_H_
#define _OPENMM_WRAPPER_H_

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef GMX_OPENMM
void* openmm_init(FILE *fplog, const char *platformOptStr,
                    t_inputrec *ir,
                    gmx_mtop_t *top_global, gmx_localtop_t *top,
                    t_mdatoms *mdatoms, t_forcerec *fr, t_state *state);

void openmm_take_one_step(void* data);

void openmm_take_steps(void* data, int nsteps);

void openmm_copy_state(void *data,
                        t_state *state, double *time,
                        rvec f[], gmx_enerdata_t *enerd,
                        gmx_bool includePos, gmx_bool includeVel, gmx_bool includeForce, gmx_bool includeEnergy);

void openmm_cleanup(FILE *fplog, void* data);
#else 
/* dummy versions of the wrapper functions to enable compilation of 
   do_md_openmm even when OpenMM is not used */ 
void* openmm_init(FILE *fplog, const char *platformOptStr,
                    t_inputrec *ir,
                    gmx_mtop_t *top_global, gmx_localtop_t *top,
                    t_mdatoms *mdatoms, t_forcerec *fr, t_state *state){return NULL;}

void openmm_take_one_step(void* data){}

void openmm_take_steps(void* data, int nsteps){}

void openmm_copy_state(void *data,
                        t_state *state, double *time,
                        rvec f[], gmx_enerdata_t *enerd,
                        gmx_bool includePos, gmx_bool includeVel, gmx_bool includeForce, gmx_bool includeEnergy){}

void openmm_cleanup(FILE *fplog, void* data){}

#endif /*GMX_OPENMM*/


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* _OPENMM_WRAPPER_H_ */

