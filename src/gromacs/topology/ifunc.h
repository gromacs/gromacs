/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_IFUNC_H
#define GMX_TOPOLOGY_IFUNC_H

#include "gromacs/fda/FDA.h"
#include "gromacs/topology/idef.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_fcdata;
struct t_graph;
struct t_mdatoms;
struct t_pbc;

/* Real vector type with an additional, unused 4th element.
 * This type is used to allow aligned 4-wide SIMD loads and stores.
 */
typedef real rvec4[4];

typedef real t_ifunc (int nbonds, const t_iatom iatoms[],
                      const t_iparams iparams[],
                      const rvec x[], rvec4 f[], rvec fshift[],
                      const struct t_pbc *pbc, const struct t_graph *g,
                      real lambda, real *dvdlambda,
                      const struct t_mdatoms *md, struct t_fcdata *fcd,
                      int *ddgatindex
#ifdef BUILD_WITH_FDA
                      , FDA *fda
#endif
                     );

/*
 * The function type t_ifunc() calculates one interaction, using iatoms[]
 * and iparams. Within the function the number of atoms to be used is
 * known. Within the function only the atomid part of the iatoms[] array
 * is supplied, not the type field (see also t_ilist). The function
 * returns the potential energy. If pbc==NULL the coordinates in x are
 * assumed to be such that no calculation of PBC is necessary,
 * If pbc!=NULL a full PBC calculation is performed.
 * If g!=NULL it is used for determining the shift forces.
 * With domain decomposition ddgatindex can be used for getting global
 * atom numbers for warnings and error messages.
 * ddgatindex is NULL when domain decomposition is not used.
 */

#define IF_NULL       0
#define IF_BOND       1
#define IF_VSITE      1<<1
#define IF_CONSTRAINT 1<<2
#define IF_CHEMBOND   1<<3
#define IF_BTYPE      1<<4
#define IF_ATYPE      1<<5
#define IF_TABULATED  1<<6
#define IF_LIMZERO    1<<7
/* These flags tell to some of the routines what can be done with this
 * item in the list.
 * With IF_BOND a bonded interaction will be calculated.
 * With IF_BTYPE grompp can convert the bond to a Morse potential.
 * With IF_BTYPE or IF_ATYPE the bond/angle can be converted to
 * a constraint or used for vsite parameter determination by grompp.
 * IF_LIMZERO indicates that for a bonded interaction the potential
 * does goes to zero for large distances, thus if such an interaction
 * it not assigned to any node by the domain decompostion, the simulation
 * still continue, if mdrun has been told so.
 */
typedef struct
{
    const char *name;         /* the name of this function			*/
    const char *longname;     /* The name for printing etc.                   */
    int         nratoms;      /* nr of atoms needed for this function		*/
    int         nrfpA, nrfpB; /* number of parameters for this function.      */
                              /* this corresponds to the number of params in  */
                              /* iparams struct! (see idef.h)                 */
    /* A and B are for normal and free energy components respectively.    */
    unsigned long   flags;    /* Flags (see above)                            */
    int             nrnb_ind; /* index for nrnb (-1 if unknown)               */
    t_ifunc        *ifunc;    /* the function it self				*/
} t_interaction_function;

#define NRFPA(ftype) (interaction_function[(ftype)].nrfpA)
#define NRFPB(ftype) (interaction_function[(ftype)].nrfpB)
#define NRFP(ftype)  (NRFPA(ftype)+NRFPB(ftype))
#define NRAL(ftype) (interaction_function[(ftype)].nratoms)

#define IS_CHEMBOND(ftype) (interaction_function[(ftype)].nratoms == 2 && (interaction_function[(ftype)].flags & IF_CHEMBOND))
/* IS_CHEMBOND tells if function type ftype represents a chemical bond */

/* IS_ANGLE tells if a function type ftype represents an angle
 * Per Larsson, 2007-11-06
 */
#define IS_ANGLE(ftype) (interaction_function[(ftype)].nratoms == 3 && (interaction_function[(ftype)].flags & IF_ATYPE))
#define IS_VSITE(ftype) (interaction_function[(ftype)].flags & IF_VSITE)

#define IS_TABULATED(ftype) (interaction_function[(ftype)].flags & IF_TABULATED)

extern const t_interaction_function interaction_function[F_NRE];
/* initialised interaction functions descriptor				*/

#ifdef __cplusplus
}
#endif

#endif
