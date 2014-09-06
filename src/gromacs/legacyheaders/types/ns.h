/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef GMX_LEGACYHEADERS_TYPES_NS_H
#define GMX_LEGACYHEADERS_TYPES_NS_H

#include "gromacs/legacyheaders/types/nblist.h"
#include "gromacs/legacyheaders/types/nsgrid.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {
    eNL_VDWQQ, eNL_VDW, eNL_QQ,
    eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE,
    eNL_VDWQQ_WATER, eNL_QQ_WATER,
    eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER,
    eNL_NR
};

#define MAX_CG 1024

typedef struct {
    int     ncg;
    int     nj;
    atom_id jcg[MAX_CG];
} t_ns_buf;

typedef struct {
    gmx_bool      bCGlist;
    atom_id      *simple_aaj;
    t_grid       *grid;
    t_excl       *bexcl;
    gmx_bool     *bHaveVdW;
    t_ns_buf    **ns_buf;
    gmx_bool     *bExcludeAlleg;
    int           nra_alloc;
    int           cg_alloc;
    atom_id     **nl_sr;
    int          *nsr;
    atom_id     **nl_lr_ljc;
    atom_id     **nl_lr_one;
    int          *nlr_ljc;
    int          *nlr_one;
    /* the nblists should probably go in here */
    gmx_bool      nblist_initialized; /* has the nblist been initialized?  */
    int           dump_nl;            /* neighbour list dump level (from env. var. GMX_DUMP_NL)*/
} gmx_ns_t;

#ifdef __cplusplus
}
#endif
#endif
