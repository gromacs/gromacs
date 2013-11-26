/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#ifndef _nblist_h
#define _nblist_h

#ifdef __cplusplus
extern "C" {
#endif


typedef unsigned long t_excl;

/* The maximum charge group size because of minimum size of t_excl
 * could be 32 bits.
 */
#define MAX_CHARGEGROUP_SIZE 32

/* The maximum charge group size for CG-CG nblists.
 * The excl entry in t_nblist uses blocks of this size.
 */
#define MAX_CGCGSIZE 32

typedef struct
{
    int             igeometry;    /* The type of list (atom, water, etc.)  */
    int             ielec;        /* Coulomb loop type index for kernels   */
    int             ielecmod;     /* Coulomb modifier (e.g. switch/shift)  */
    int             ivdw;         /* VdW loop type index for kernels       */
    int             ivdwmod;      /* VdW modifier (e.g. switch/shift)      */
    int             type;         /* Type of interaction, listed in
                                     gmx_nblist_interaction_type           */

    int             nri, maxnri;  /* Current/max number of i particles	   */
    int             nrj, maxnrj;  /* Current/max number of j particles	   */
    int *           iinr;         /* The i-elements                        */
    int *           iinr_end;     /* The end atom, only with enlistCG      */
    int *           gid;          /* Index in energy arrays                */
    int *           shift;        /* Shift vector index                    */
    int *           jindex;       /* Index in jjnr                         */
    int *           jjnr;         /* The j-atom list                       */
    int *           jjnr_end;     /* The end atom, only with enltypeCG     */
    char *          excl_fep;     /* Exclusions for FEP with Verlet scheme */
    t_excl *        excl;         /* Exclusions, only with enltypeCG       */

    /* We use separate pointers for kernels that compute both potential
     * and force (vf suffix), only potential (v) or only force (f)
     */
    void *          kernelptr_vf;
    void *          kernelptr_v;
    void *          kernelptr_f;

    /* Pad the list of neighbors for each i atom with "-1" entries up to the
     * simd_padding_width, if it is larger than 0. This is necessary for many
     * accelerated kernels using single-instruction multiple-data operations
     * internally.
     */
    int             simd_padding_width;

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

#ifdef __cplusplus
}
#endif

#endif
