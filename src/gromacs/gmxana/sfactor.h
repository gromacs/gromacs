/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_GMXANA_SFACTOR_H
#define GMX_GMXANA_SFACTOR_H

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;
struct t_topology;
struct t_trxframe;

typedef struct gmx_structurefactors gmx_structurefactors_t;

typedef struct structure_factor structure_factor_t;

typedef struct reduced_atom reduced_atom_t;

int* create_indexed_atom_type(reduced_atom_t* atm, int size);

void compute_structure_factor(structure_factor_t* sft,
                              matrix              box,
                              reduced_atom_t*     red,
                              int                 isize,
                              real                start_q,
                              real                end_q,
                              int                 group,
                              real**              sf_table);

gmx_structurefactors_t* gmx_structurefactors_init(const char* datfn);

void gmx_structurefactors_done(gmx_structurefactors_t* gsf);

int gmx_structurefactors_get_sf(gmx_structurefactors_t* gsf, int elem, real a[4], real b[4], real* c);

real** gmx_structurefactors_table(gmx_structurefactors_t* gsf, real momentum, real ref_k, real lambda, int n_angles);

void save_data(structure_factor_t* sft, const char* file, int ngrps, real start_q, real end_q, const gmx_output_env_t* oenv);

double CMSF(gmx_structurefactors_t* gsf, int type, int nh, double lambda, double sin_theta);

int return_atom_type(const char* name, gmx_structurefactors_t* gsf);

void rearrange_atoms(reduced_atom_t*         positions,
                     struct t_trxframe*      fr,
                     const int*              index,
                     int                     isize,
                     const t_topology*       top,
                     bool                    flag,
                     gmx_structurefactors_t* gsf);

int do_scattering_intensity(const char*             fnTPS,
                            const char*             fnNDX,
                            const char*             fnXVG,
                            const char*             fnTRX,
                            const char*             fnDAT,
                            real                    start_q,
                            real                    end_q,
                            real                    energy,
                            int                     ng,
                            const gmx_output_env_t* oenv);

t_complex*** rc_tensor_allocation(int x, int y, int z);

real** compute_scattering_factor_table(gmx_structurefactors_t* gsf, structure_factor_t* sft);

#endif
