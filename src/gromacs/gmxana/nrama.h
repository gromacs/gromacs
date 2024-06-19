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

#ifndef GMX_GMXANA_NRAMA_H
#define GMX_GMXANA_NRAMA_H

#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;
enum class PbcType : int;
struct t_idef;

typedef struct
{
    gmx_bool bShow;
    char*    label;
    int      iphi, ipsi; /* point in the dih array of xr... */
} t_phipsi;

typedef struct
{
    int  ai[4];
    int  mult;
    real phi0;
    real ang;
} t_dih;

typedef struct
{
    int               ndih;
    t_dih*            dih;
    int               npp;
    t_phipsi*         pp;
    t_trxstatus*      traj;
    int               natoms;
    int               amin, amax;
    real              t;
    rvec*             x;
    matrix            box;
    t_idef*           idef;
    PbcType           pbcType;
    gmx_output_env_t* oenv;
} t_xrama;

t_topology* init_rama(gmx_output_env_t* oenv, const char* infile, const char* topfile, t_xrama* xr, int mult);

gmx_bool new_data(t_xrama* xr);

#endif /* GMX_GMXANA_NRAMA_H */
