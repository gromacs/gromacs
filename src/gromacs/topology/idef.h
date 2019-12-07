/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_IDEF_H
#define GMX_TOPOLOGY_IDEF_H

#include <cstdio>

#include <array>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

typedef union t_iparams {
    /* Some parameters have A and B values for free energy calculations.
     * The B values are not used for regular simulations of course.
     * Free Energy for nonbondeds can be computed by changing the atom type.
     * The harmonic type is used for all harmonic potentials:
     * bonds, angles and improper dihedrals
     */
    struct
    {
        real a, b, c;
    } bham;
    struct
    {
        real rA, krA, rB, krB;
    } harmonic;
    struct
    {
        real klinA, aA, klinB, aB;
    } linangle;
    struct
    {
        real lowA, up1A, up2A, kA, lowB, up1B, up2B, kB;
    } restraint;
    /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */
    struct
    {
        real b0, kb, kcub;
    } cubic;
    struct
    {
        real bm, kb;
    } fene;
    struct
    {
        real r1e, r2e, krr;
    } cross_bb;
    struct
    {
        real r1e, r2e, r3e, krt;
    } cross_ba;
    struct
    {
        real thetaA, kthetaA, r13A, kUBA, thetaB, kthetaB, r13B, kUBB;
    } u_b;
    struct
    {
        real theta, c[5];
    } qangle;
    struct
    {
        real alpha;
    } polarize;
    struct
    {
        real alpha, drcut, khyp;
    } anharm_polarize;
    struct
    {
        real al_x, al_y, al_z, rOH, rHH, rOD;
    } wpol;
    struct
    {
        real a, alpha1, alpha2, rfac;
    } thole;
    struct
    {
        real c6, c12;
    } lj;
    struct
    {
        real c6A, c12A, c6B, c12B;
    } lj14;
    struct
    {
        real fqq, qi, qj, c6, c12;
    } ljc14;
    struct
    {
        real qi, qj, c6, c12;
    } ljcnb;
    /* Proper dihedrals can not have different multiplicity when
     * doing free energy calculations, because the potential would not
     * be periodic anymore.
     */
    struct
    {
        real phiA, cpA;
        int  mult;
        real phiB, cpB;
    } pdihs;
    struct
    {
        real dA, dB;
    } constr;
    /* Settle can not be used for Free energy calculations of water bond geometry.
     * Use shake (or lincs) instead if you have to change the water bonds.
     */
    struct
    {
        real doh, dhh;
    } settle;
    struct
    {
        real b0A, cbA, betaA, b0B, cbB, betaB;
    } morse;
    struct
    {
        real pos0A[DIM], fcA[DIM], pos0B[DIM], fcB[DIM];
    } posres;
    struct
    {
        real pos0[DIM], r, k;
        int  geom;
    } fbposres;
    struct
    {
        real rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];
    } rbdihs;
    struct
    {
        real cbtcA[NR_CBTDIHS], cbtcB[NR_CBTDIHS];
    } cbtdihs;
    struct
    {
        real a, b, c, d, e, f;
    } vsite;
    struct
    {
        int  n;
        real a;
    } vsiten;
    /* NOTE: npair is only set after reading the tpx file */
    struct
    {
        real low, up1, up2, kfac;
        int  type, label, npair;
    } disres;
    struct
    {
        real phiA, dphiA, kfacA, phiB, dphiB, kfacB;
    } dihres;
    struct
    {
        int  ex, power, label;
        real c, obs, kfac;
    } orires;
    struct
    {
        int  table;
        real kA;
        real kB;
    } tab;
    struct
    {
        int cmapA, cmapB;
    } cmap;
    struct
    {
        real buf[MAXFORCEPARAM];
    } generic; /* Conversion */
} t_iparams;

typedef int t_functype;

/* List of listed interactions, see description further down.
 *
 * TODO: Consider storing the function type as well.
 * TODO: Consider providing per interaction access.
 */
struct InteractionList
{
    /* Returns the total number of elements in iatoms */
    int size() const { return gmx::ssize(iatoms); }

    /* List of interactions, see explanation further down */
    std::vector<int> iatoms;
};

/* List of interaction lists, one list for each interaction type
 *
 * TODO: Consider only including entries in use instead of all F_NRE
 */
typedef std::array<InteractionList, F_NRE> InteractionLists;

/* Deprecated list of listed interactions.
 *
 * The nonperturbed/perturbed interactions are now separated (sorted) in the
 * ilist, such that the first 0..(nr_nonperturbed-1) ones are exactly that, and
 * the remaining ones from nr_nonperturbed..(nr-1) are perturbed bonded
 * interactions.
 */
struct t_ilist
{
    /* Returns the total number of elements in iatoms */
    int size() const { return nr; }

    int      nr;
    int      nr_nonperturbed;
    t_iatom* iatoms;
    int      nalloc;
};

/* TODO: Replace t_ilist in gmx_localtop_t by InteractionList.
 *       The nr_nonperturbed functionality needs to be ported.
 *       Remove t_topology.
 *       Remove t_ilist and remove templating on list type
 *       in mshift.cpp, constr.cpp, vsite.cpp and domdec_topology.cpp.
 */

/*
 * The structs InteractionList and t_ilist defines a list of atoms with their interactions.
 * General field description:
 *   int nr
 *      the size (nr elements) of the interactions array (iatoms[]).
 *   t_iatom *iatoms
 *  specifies which atoms are involved in an interaction of a certain
 *       type. The layout of this array is as follows:
 *
 *        +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
 *        |type1|at1|at2|at3|type2|at1|at2|type1|at1|at2|at3|type3|at1|at2|
 *        +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
 *
 *  So for interaction type type1 3 atoms are needed, and for type2 and
 *      type3 only 2. The type identifier is used to select the function to
 *      calculate the interaction and its actual parameters. This type
 *      identifier is an index in a params[] and functype[] array.
 */

/*! \brief Type for returning a list of InteractionList references
 *
 * TODO: Remove when the function type is made part of InteractionList
 */
struct InteractionListHandle
{
    const int               functionType; //!< The function type
    const std::vector<int>& iatoms;       //!< Reference to interaction list
};

/*! \brief Returns a list of all non-empty InteractionList entries with any of the interaction flags in \p flags set
 *
 * \param[in] ilists  Set of interaction lists
 * \param[in] flags   Bit mask with one or more IF_... bits set
 */
static inline std::vector<InteractionListHandle> extractILists(const InteractionLists& ilists, int flags)
{
    std::vector<InteractionListHandle> handles;
    for (size_t ftype = 0; ftype < ilists.size(); ftype++)
    {
        if ((interaction_function[ftype].flags & flags) && ilists[ftype].size() > 0)
        {
            handles.push_back({ static_cast<int>(ftype), ilists[ftype].iatoms });
        }
    }
    return handles;
}

/*! \brief Returns the stride for the iatoms array in \p ilistHandle
 *
 * \param[in] ilistHandle  The ilist to return the stride for
 */
static inline int ilistStride(const InteractionListHandle& ilistHandle)
{
    return 1 + NRAL(ilistHandle.functionType);
}

struct gmx_cmapdata_t
{
    std::vector<real> cmap; /* Has length 4*grid_spacing*grid_spacing, */
    /* there are 4 entries for each cmap type (V,dVdx,dVdy,d2dVdxdy) */
};

struct gmx_cmap_t
{
    int                         grid_spacing = 0; /* Grid spacing */
    std::vector<gmx_cmapdata_t> cmapdata; /* Lists of grids with actual, pre-interpolated data */
};


enum
{
    ilsortUNKNOWN,
    ilsortNO_FE,
    ilsortFE_UNSORTED,
    ilsortFE_SORTED
};

typedef struct t_idef
{
    int         ntypes;
    int         atnr;
    t_functype* functype;
    t_iparams*  iparams;
    real        fudgeQQ;
    gmx_cmap_t* cmap_grid;
    t_iparams * iparams_posres, *iparams_fbposres;
    int         iparams_posres_nalloc, iparams_fbposres_nalloc;

    t_ilist il[F_NRE];
    int     ilsort;
} t_idef;

/*
 * The struct t_idef defines all the interactions for the complete
 * simulation. The structure is setup in such a way that the multinode
 * version of the program  can use it as easy as the single node version.
 * General field description:
 *   int ntypes
 *      defines the number of elements in functype[] and param[].
 *   int nodeid
 *      the node id (if parallel machines)
 *   int atnr
 *      the number of atomtypes
 *   t_functype *functype
 *      array of length ntypes, defines for every force type what type of
 *      function to use. Every "bond" with the same function but different
 *      force parameters is a different force type. The type identifier in the
 *      forceatoms[] array is an index in this array.
 *   t_iparams *iparams
 *      array of length ntypes, defines the parameters for every interaction
 *      type. The type identifier in the actual interaction list
 *      (ilist[ftype].iatoms[]) is an index in this array.
 *   gmx_cmap_t cmap_grid
 *      the grid for the dihedral pair correction maps.
 *   t_iparams *iparams_posres, *iparams_fbposres
 *      defines the parameters for position restraints only.
 *      Position restraints are the only interactions that have different
 *      parameters (reference positions) for different molecules
 *      of the same type. ilist[F_POSRES].iatoms[] is an index in this array.
 *   t_ilist il[F_NRE]
 *      The list of interactions for each type. Note that some,
 *      such as LJ and COUL will have 0 entries.
 *   int ilsort
 *      The state of the sorting of il, values are provided above.
 */

namespace gmx
{
class TextWriter;
} // namespace gmx

void printInteractionParameters(gmx::TextWriter* writer, t_functype ftype, const t_iparams& iparams);
void pr_iparams(FILE* fp, t_functype ftype, const t_iparams& iparams);
void pr_ilist(FILE*                  fp,
              int                    indent,
              const char*            title,
              const t_functype*      functype,
              const InteractionList& ilist,
              gmx_bool               bShowNumbers,
              gmx_bool               bShowParameters,
              const t_iparams*       iparams);
void pr_idef(FILE* fp, int indent, const char* title, const t_idef* idef, gmx_bool bShowNumbers, gmx_bool bShowParameters);

/*! \brief
 * Properly initialize idef struct.
 *
 * \param[in] idef Pointer to idef struct to initialize.
 */
void init_idef(t_idef* idef);

/*! \brief
 * Properly clean up idef struct.
 *
 * \param[in] idef Pointer to idef struct to clean up.
 */
void done_idef(t_idef* idef);

void copy_ilist(const t_ilist* src, t_ilist* dst);

#endif
