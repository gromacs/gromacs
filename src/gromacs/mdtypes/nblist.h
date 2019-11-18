/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2019, by the GROMACS development team, led by
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
#ifndef GMX_MDTYPES_NBLIST_H
#define GMX_MDTYPES_NBLIST_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

typedef unsigned long t_excl;

/* The interactions contained in a (possibly merged) table
 * for computing electrostatic, VDW repulsion and/or VDW dispersion
 * contributions.
 */
enum gmx_table_interaction
{
    GMX_TABLE_INTERACTION_ELEC,
    GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWDISP,
    GMX_TABLE_INTERACTION_NR
};

/* Different formats for table data. Cubic spline tables are typically stored
 * with the four Y,F,G,H intermediate values (check tables.c for format), which
 * makes it easy to load with a single 4-way SIMD instruction too.
 * Linear tables only need one value per table point, or two if both V and F
 * are calculated. However, with SIMD instructions this makes the loads unaligned,
 * and in that case we store the data as F, D=F(i+1)-F(i), V, and then a blank value,
 * which again makes it possible to load as a single instruction.
 */
enum gmx_table_format
{
    GMX_TABLE_FORMAT_CUBICSPLINE_YFGH,
    GMX_TABLE_FORMAT_LINEAR_VF,
    GMX_TABLE_FORMAT_LINEAR_V,
    GMX_TABLE_FORMAT_LINEAR_F,
    GMX_TABLE_FORMAT_LINEAR_FDV0,
    GMX_TABLE_FORMAT_NR
};

enum
{
    eNL_VDWQQ,
    eNL_VDW,
    eNL_QQ,
    eNL_VDWQQ_FREE,
    eNL_VDW_FREE,
    eNL_QQ_FREE,
    eNL_VDWQQ_WATER,
    eNL_QQ_WATER,
    eNL_VDWQQ_WATERWATER,
    eNL_QQ_WATERWATER,
    eNL_NR
};

#define MAX_CG 1024

typedef struct
{
    int ncg;
    int nj;
    int jcg[MAX_CG];
} t_ns_buf;


/* The maximum charge group size because of minimum size of t_excl
 * could be 32 bits.
 */
#define MAX_CHARGEGROUP_SIZE 32

/* The maximum charge group size for CG-CG nblists.
 * The excl entry in t_nblist uses blocks of this size.
 */
#define MAX_CGCGSIZE 32

typedef struct t_nblist
{
    int igeometry; /* The type of list (atom, water, etc.)  */
    int ielec;     /* Coulomb loop type index for kernels   */
    int ielecmod;  /* Coulomb modifier (e.g. switch/shift)  */
    int ivdw;      /* VdW loop type index for kernels       */
    int ivdwmod;   /* VdW modifier (e.g. switch/shift)      */
    int type;      /* Type of interaction, listed in
                      gmx_nblist_interaction_type           */

    int     nri, maxnri; /* Current/max number of i particles	   */
    int     nrj, maxnrj; /* Current/max number of j particles	   */
    int*    iinr;        /* The i-elements                        */
    int*    iinr_end;    /* The end atom, only with enlistCG      */
    int*    gid;         /* Index in energy arrays                */
    int*    shift;       /* Shift vector index                    */
    int*    jindex;      /* Index in jjnr                         */
    int*    jjnr;        /* The j-atom list                       */
    int*    jjnr_end;    /* The end atom, only with enltypeCG     */
    char*   excl_fep;    /* Exclusions for FEP with Verlet scheme */
    t_excl* excl;        /* Exclusions, only with enltypeCG       */

    /* We use separate pointers for kernels that compute both potential
     * and force (vf suffix), only potential (v) or only force (f)
     */
    void* kernelptr_vf;
    void* kernelptr_v;
    void* kernelptr_f;

    /* Pad the list of neighbors for each i atom with "-1" entries up to the
     * simd_padding_width, if it is larger than 0. This is necessary for many
     * accelerated kernels using single-instruction multiple-data operations
     * internally.
     */
    int simd_padding_width;

} t_nblist;

/* For atom I =  nblist->iinr[N] (0 <= N < nblist->nri) there can be
 * several neighborlists (N's), for different energy groups (gid) and
 * different shifts (shift).
 * For corresponding J atoms for each list start at:
 * nblist->jjnr[JI]
 * with nblist->jindex[N] <= JI < nblist->jindex[N+1]
 *
 * enlist is of the form enlistUNIT1_UNIT2:
 * UNIT ATOM:  there is one atom: iinr[N] or jjnr[JI]
 * UNIT SPC:   there are 3 atoms: iinr[N],iinr[N]+1,iinr[N]+2, jjnr analog.
 * UNIT TIP4P: there are 4 atoms: iinr[N],...,iinr[N]+3, jjnr analog.
 * UNIT CG:    there are N atoms: iinr[N],...,iinr_end[N]-1, jjnr analog.
 *
 * Clear?
 */

/* Structure describing the data in a single table */
struct t_forcetable
{
    t_forcetable(enum gmx_table_interaction interaction, enum gmx_table_format format);

    ~t_forcetable();

    enum gmx_table_interaction interaction; /* Types of interactions stored in this table */
    enum gmx_table_format      format;      /* Interpolation type and data format */

    real  r;     /* range of the table */
    int   n;     /* n+1 is the number of table points */
    real  scale; /* distance (nm) between two table points */
    real* data;  /* the actual table data */

    /* Some information about the table layout. This can also be derived from the interpolation
     * type and the table interactions, but it is convenient to have here for sanity checks, and it
     * makes it much easier to access the tables in the nonbonded kernels when we can set the data
     * from variables. It is always true that stride = formatsize*ninteractions
     */
    int formatsize; /* Number of fp variables for each table point (1 for F, 2 for VF, 4 for YFGH, etc.) */
    int ninteractions; /* Number of interactions in table, 1 for coul-only, 3 for coul+rep+disp. */
    int stride; /* Distance to next table point (number of fp variables per table point in total) */
};

struct t_nblists
{
    struct t_forcetable* table_elec;
    struct t_forcetable* table_vdw;
    struct t_forcetable* table_elec_vdw;

    /* The actual neighbor lists, short and long range, see enum above
     * for definition of neighborlist indices.
     */
    struct t_nblist nlist_sr[eNL_NR];
    struct t_nblist nlist_lr[eNL_NR];
};

struct gmx_ns_t
{
    gmx_bool       bCGlist;
    int*           simple_aaj;
    struct t_grid* grid;
    t_excl*        bexcl;
    gmx_bool*      bHaveVdW;
    t_ns_buf**     ns_buf;
    gmx_bool*      bExcludeAlleg;
    int            nra_alloc;
    int            cg_alloc;
    int**          nl_sr;
    int*           nsr;
    int**          nl_lr_ljc;
    int**          nl_lr_one;
    int*           nlr_ljc;
    int*           nlr_one;
    /* the nblists should probably go in here */
    gmx_bool nblist_initialized; /* has the nblist been initialized?  */
    int      dump_nl;            /* neighbour list dump level (from env. var. GMX_DUMP_NL)*/
};

#endif /* GMX_MDTYPES_NBLIST_H */
