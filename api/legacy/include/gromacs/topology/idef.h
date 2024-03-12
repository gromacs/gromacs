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
/*! \file
 * \brief
 * Defines interaction types and parameters.
 *
 * \ingroup module_topology
 * \inlibraryapi
 */
#ifndef GMX_TOPOLOGY_IDEF_H
#define GMX_TOPOLOGY_IDEF_H

#include <cstdio>

#include <array>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/real.h"

struct gmx_ffparams_t;

/*! \brief
 * Union of various interaction paramters.
 *
 * Some parameters have A and B values for free energy calculations.
 * The B values are not used for regular simulations of course.
 * Free Energy for nonbondeds can be computed by changing the atom type.
 * The harmonic type is used for all harmonic potentials:
 * bonds, angles and improper dihedrals.
 *
 * No free energy supported for cubic bonds, FENE, WPOL or cross terms.
 */
typedef union t_iparams
{
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
        real a, alpha1, alpha2;
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
    /*! \brief Proper dihedrals.
     *
     * Proper dihedrals can not have different multiplicity when
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
    /*! \brief SETTLE.
     *
     * SETTLE can not be used for Free energy calculations of water bond geometry.
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

/*! \brief List of listed interactions.
 *
 * Defines a list of atoms with their interactions.
 *
 * TODO: Consider storing the function type as well.
 * TODO: Consider providing per interaction access.
 */
struct InteractionList
{
    //! \return The total number of elements in iatoms.
    int size() const { return static_cast<int>(iatoms.size()); }

    //! \return Whether the list is empty.
    bool empty() const { return iatoms.empty(); }

    /*! \brief Adds one interaction to the list.
     *
     * \param[in] parameterType Type.
     * \param[in] atoms Array of atoms indices.
     */
    template<std::size_t numAtoms>
    void push_back(const int parameterType, const std::array<int, numAtoms>& atoms)
    {
        const std::size_t oldSize = iatoms.size();
        iatoms.resize(iatoms.size() + 1 + numAtoms);
        iatoms[oldSize] = parameterType;
        for (std::size_t i = 0; i < numAtoms; i++)
        {
            iatoms[oldSize + 1 + i] = atoms[i];
        }
    }

    /*! \brief Adds one interaction to the list.
     *
     * \param[in] parameterType Type.
     * \param[in] numAtoms Number of atoms.
     * \param[in] atoms Atoms indices.
     */
    void push_back(const int parameterType, const int numAtoms, const int* atoms)
    {
        const std::size_t oldSize = iatoms.size();
        iatoms.resize(iatoms.size() + 1 + numAtoms);
        iatoms[oldSize] = parameterType;
        for (int i = 0; i < numAtoms; i++)
        {
            iatoms[oldSize + 1 + i] = atoms[i];
        }
    }

    /*! \brief Appends interaction list \p ilist at the back of the list.
     *
     * \param[in] ilist Interaction list to append.
     */
    void append(const InteractionList& ilist)
    {
        iatoms.insert(iatoms.end(), ilist.iatoms.begin(), ilist.iatoms.end());
    }

    /*! \brief Clears the list.
     */
    void clear() { iatoms.clear(); }

    /*! \brief List of interactions.
     *
     * Specifies which atoms are involved in an interaction of a certain type.
     * The layout of this array is as follows:
     *
     *     +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
     *     |type1|at1|at2|at3|type2|at1|at2|type1|at1|at2|at3|type3|at1|at2|
     *     +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
     *
     * So for interaction type type1 3 atoms are needed, and for type2 and
     * type3 only 2. The type identifier is used to select the function to
     * calculate the interaction and its actual parameters. This type
     * identifier is an index in a `params[]` and `functype[]` array.
     */
    std::vector<int> iatoms;
};

/*! \brief List of interaction lists.
 *
 * Contains one list for each interaction type.
 *
 * TODO: Consider only including entries in use instead of all F_NRE
 */
using InteractionLists = std::array<InteractionList, F_NRE>;

/*! \brief Deprecated list of listed interactions.
 */
struct t_ilist
{
    //! \return The total number of elements in iatoms.
    int size() const { return nr; }

    //! \return Whether the list is empty.
    bool empty() const { return nr == 0; }

    //! The size (nr elements) of the interactions array (iatoms[]).
    int      nr;
    t_iatom* iatoms;
    int      nalloc;
};

/* TODO: Remove t_ilist and remove templating on list type in mshift.cpp */

/*! \brief Type for returning a list of InteractionList references.
 *
 * TODO: Remove when the function type is made part of InteractionList.
 */
struct InteractionListHandle
{
    const int               functionType; //!< The function type.
    const std::vector<int>& iatoms;       //!< Reference to interaction list.
};

/*! \brief Returns a list of all non-empty InteractionList entries with any of the interaction flags in \p flags set
 *
 * \param[in] ilists Set of interaction lists.
 * \param[in] flags Bit mask with one or more IF_... bits set.
 */
static inline std::vector<InteractionListHandle> extractILists(const InteractionLists& ilists, int flags)
{
    std::vector<InteractionListHandle> handles;
    for (size_t ftype = 0; ftype < ilists.size(); ftype++)
    {
        if ((interaction_function[ftype].flags & flags) && !ilists[ftype].empty())
        {
            handles.push_back({ static_cast<int>(ftype), ilists[ftype].iatoms });
        }
    }
    return handles;
}

/*! \brief Returns the stride for the iatoms array in \p ilistHandle.
 *
 * \param[in] ilistHandle The ilist to return the stride for.
 */
static inline int ilistStride(const InteractionListHandle& ilistHandle)
{
    return 1 + NRAL(ilistHandle.functionType);
}

/*! \brief CMAP data for a single grid. */
struct gmx_cmapdata_t
{
    /*! Has length \f$4 \times grid\_spacing \times grid\_spacing\f$,
     * there are 4 entries for each CMAP type \f$(V, dV dx, dV dy, d^2 dV dx dy)\f$
     */
    std::vector<real> cmap;
};

/*! \brief CMAP grids. */
struct gmx_cmap_t
{
    int                         grid_spacing = 0; //!< Grid spacing.
    std::vector<gmx_cmapdata_t> cmapdata; //!< Lists of grids with actual, pre-interpolated data.
};

/*! \brief Sorting of the interaction list. */
enum
{
    ilsortUNKNOWN,  //!< Sorting unknown.
    ilsortNO_FE,    //!< Sorting without the free energy interactions.
    ilsortFE_SORTED //!< Sorting with the free energy interactions also sorted.
};

/*! \brief
 * Struct with list of interaction parameters and lists of interactions
 *
 * TODO: Convert to a proper class with private data members so we can
 * ensure that the free-energy sorting and sorting setting is consistent.
 */
class InteractionDefinitions
{
public:
    /*! \brief Constructor
     *
     * \param[in] ffparams The interaction parameters, the lifetime of the created object should not exceed the lifetime of the passed parameters.
     */
    InteractionDefinitions(const gmx_ffparams_t& ffparams);

    //! \brief Clears data not read in from ffparams.
    void clear();

    //! The interaction parameters.
    const std::vector<t_iparams>& iparams;
    //! The function type per type.
    const std::vector<int>& functype;
    //! Position restraint interaction parameters.
    std::vector<t_iparams> iparams_posres;
    //! Flat-bottomed position restraint parameters.
    std::vector<t_iparams> iparams_fbposres;
    //! The list of interactions for each type. Note that some, such as LJ and COUL will have 0 entries.
    std::array<InteractionList, F_NRE> il;
    //! The number of non-perturbed interactions at the start of each entry in il.
    std::array<int, F_NRE> numNonperturbedInteractions;
    //! The sorting state of interaction in il.
    int ilsort = ilsortUNKNOWN;
    //! The dihedral correction maps.
    gmx_cmap_t cmap_grid;
};

/*! \brief Deprecated interation definitions, used in t_topology
 *
 * The struct t_idef defines all the interactions for the complete
 * simulation. The structure is setup in such a way that the multinode
 * version of the program  can use it as easy as the single node version.
 */
struct t_idef
{
    //! Defines the number of elements in functype[] and param[].
    int ntypes;
    //! The number of atomtypes.
    int atnr;
    /*! Array of length ntypes, defines for every force type what type of
     * function to use. Every "bond" with the same function but different
     * force parameters is a different force type. The type identifier in the
     * forceatoms[] array is an index in this array.
     */
    t_functype* functype;
    /*! Array of length ntypes, defines the parameters for every interaction
     * type. The type identifier in the actual interaction list
     * (ilist[ftype].iatoms[]) is an index in this array.
     */
    t_iparams* iparams;
    real       fudgeQQ;
    /*! Defines the parameters for position restraints only.
     * Position restraints are the only interactions that have different
     * parameters (reference positions) for different molecules
     * of the same type. ilist[F_POSRES].iatoms[] is an index in this array.
     */
    t_iparams *iparams_posres, *iparams_fbposres;
    /*! The list of interactions for each type. Note that some,
     * such as LJ and COUL will have 0 entries.
     */
    t_ilist il[F_NRE];
    //! The state of the sorting of il, values are provided above.
    int ilsort;
};

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
              bool                   bShowNumbers,
              bool                   bShowParameters,
              const t_iparams*       iparams);
void pr_idef(FILE* fp, int indent, const char* title, const t_idef* idef, bool bShowNumbers, bool bShowParameters);

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

#endif
