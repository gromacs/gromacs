/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
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

