/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_TOPOLOGY_H
#define GMX_TOPOLOGY_TOPOLOGY_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/unique_cptr.h"

enum class SimulationAtomGroupType : int
{
    TemperatureCoupling,
    EnergyOutput,
    Acceleration,
    Freeze,
    User1,
    User2,
    MassCenterVelocityRemoval,
    CompressedPositionOutput,
    OrientationRestraintsFit,
    QuantumMechanics,
    Count
};

//! Short strings used for describing atom groups in log and energy files
const char* shortName(SimulationAtomGroupType type);

// const char *shortName(int type); // if necessary

/*! \brief Molecules type data: atoms, interactions and exclusions */
struct gmx_moltype_t
{
    gmx_moltype_t();

    ~gmx_moltype_t();

    /*! \brief Deleted copy assignment operator to avoid (not) freeing pointers */
    gmx_moltype_t& operator=(const gmx_moltype_t&) = delete;

    /*! \brief Default copy constructor */
    gmx_moltype_t(const gmx_moltype_t&) = default;

    char**           name;  /**< Name of the molecule type            */
    t_atoms          atoms; /**< The atoms in this molecule           */
    InteractionLists ilist; /**< Interaction list with local indices  */
    t_blocka         excls; /**< The exclusions                       */
};

/*! \brief Block of molecules of the same type, used in gmx_mtop_t */
struct gmx_molblock_t
{
    int                    type = -1; /**< The molecule type index in mtop.moltype  */
    int                    nmol = 0;  /**< The number of molecules in this block    */
    std::vector<gmx::RVec> posres_xA; /**< Position restraint coordinates for top A */
    std::vector<gmx::RVec> posres_xB; /**< Position restraint coordinates for top B */
};

/*! \brief Indices for a gmx_molblock_t, derived from other gmx_mtop_t contents */
struct MoleculeBlockIndices
{
    int numAtomsPerMolecule; /**< Number of atoms in a molecule in the block */
    int globalAtomStart;     /**< Global atom index of the first atom in the block */
    int globalAtomEnd;       /**< Global atom index + 1 of the last atom in the block */
    int globalResidueStart;  /**< Global residue index of the first residue in the block */
    int residueNumberStart; /**< Residue numbers start from this value if the number of residues per molecule is <= maxres_renum */
    int moleculeIndexStart; /**< Global molecule indexing starts from this value */
};

/*! \brief Contains the simulation atom groups.
 *
 * Organized as containers for the different
 * SimulationAtomGroupTypes
 */
struct SimulationGroups
{
    //! Groups of particles
    gmx::EnumerationArray<SimulationAtomGroupType, AtomGroupIndices> groups;
    //! Names of groups, stored as pointer to the entries in the symbol table.
    std::vector<char**> groupNames;
    //! Group numbers for the different SimulationAtomGroupType groups.
    gmx::EnumerationArray<SimulationAtomGroupType, std::vector<unsigned char>> groupNumbers;

    /*! \brief
     * Number of group numbers for a single SimulationGroup.
     *
     * \param[in] group Integer value for the group type.
     */
    int numberOfGroupNumbers(SimulationAtomGroupType group) const
    {
        return gmx::ssize(groupNumbers[group]);
    }
};

/*! \brief
 * Returns group number of an input group for a given atom.
 *
 * Returns the group \p type for \p atom in \p group, or 0 if the
 * entries for all atoms in the group are 0 and the pointer is thus null.
 *
 * \param[in] group Group to check.
 * \param[in] type  Type of group to check.
 * \param[in] atom  Atom to check if it has an entry.
 */
int getGroupType(const SimulationGroups& group, SimulationAtomGroupType type, int atom);

/* The global, complete system topology struct, based on molecule types.
 * This structure should contain no data that is O(natoms) in memory.
 *
 * TODO: Find a solution for ensuring that the derived data is in sync
 *       with the primary data, possibly by converting to a class.
 */
struct gmx_mtop_t //NOLINT(clang-analyzer-optin.performance.Padding)
{
    gmx_mtop_t();

    ~gmx_mtop_t();

    //! Name of the topology.
    char** name = nullptr;
    //! Force field parameters used.
    gmx_ffparams_t ffparams;
    //! Vector of different molecule types.
    std::vector<gmx_moltype_t> moltype;
    //! Vector of different molecule blocks.
    std::vector<gmx_molblock_t> molblock;
    //! Are there intermolecular interactions?
    bool bIntermolecularInteractions = false;
    /* \brief
     * List of intermolecular interactions using system wide
     * atom indices, either NULL or size F_NRE
     */
    std::unique_ptr<InteractionLists> intermolecular_ilist = nullptr;
    //! Number of global atoms.
    int natoms = 0;
    //! Parameter for residue numbering.
    int maxres_renum = 0;
    //! The maximum residue number in moltype
    int maxresnr = -1;
    //! Atomtype properties
    t_atomtypes atomtypes;
    //! Groups of atoms for different purposes
    SimulationGroups groups;
    //! The symbol table
    t_symtab symtab;
    //! Tells whether we have valid molecule indices
    bool haveMoleculeIndices = false;
    /*! \brief List of global atom indices of atoms between which
     * non-bonded interactions must be excluded.
     */
    std::vector<int> intermolecularExclusionGroup;

    /* Derived data  below */
    //! Indices for each molblock entry for fast lookup of atom properties
    std::vector<MoleculeBlockIndices> moleculeBlockIndices;
};

/*! \brief
 * The fully written out topology for a domain over its lifetime
 *
 * Also used in some analysis code.
 */
struct gmx_localtop_t
{
    //! Constructor used for normal operation, manages own resources.
    gmx_localtop_t();

    ~gmx_localtop_t();

    //! The interaction function definition
    t_idef idef;
    //! Atomtype properties
    t_atomtypes atomtypes;
    //! The exclusions
    t_blocka excls;
    //! Flag for domain decomposition so we don't free already freed memory.
    bool useInDomainDecomp_ = false;
};

/* The old topology struct, completely written out, used in analysis tools */
typedef struct t_topology
{
    char**      name;                        /* Name of the topology                 */
    t_idef      idef;                        /* The interaction function definition  */
    t_atoms     atoms;                       /* The atoms                            */
    t_atomtypes atomtypes;                   /* Atomtype properties                  */
    t_block     mols;                        /* The molecules                        */
    gmx_bool    bIntermolecularInteractions; /* Inter.mol. int. ?   */
    t_blocka    excls;                       /* The exclusions                       */
    t_symtab    symtab;                      /* The symbol table                     */
} t_topology;

void init_top(t_topology* top);
void done_top(t_topology* top);
// Frees both t_topology and gmx_mtop_t when the former has been created from
// the latter.
void done_top_mtop(t_topology* top, gmx_mtop_t* mtop);

bool gmx_mtop_has_masses(const gmx_mtop_t* mtop);
bool gmx_mtop_has_charges(const gmx_mtop_t* mtop);
bool gmx_mtop_has_perturbed_charges(const gmx_mtop_t& mtop);
bool gmx_mtop_has_atomtypes(const gmx_mtop_t* mtop);
bool gmx_mtop_has_pdbinfo(const gmx_mtop_t* mtop);

void pr_mtop(FILE* fp, int indent, const char* title, const gmx_mtop_t* mtop, gmx_bool bShowNumbers, gmx_bool bShowParameters);
void pr_top(FILE* fp, int indent, const char* title, const t_topology* top, gmx_bool bShowNumbers, gmx_bool bShowParameters);

/*! \brief Compare two mtop topologies.
 *
 * \param[in] fp File pointer to write to.
 * \param[in] mtop1 First topology to compare.
 * \param[in] mtop2 Second topology to compare.
 * \param[in] relativeTolerance Relative tolerance for comparison.
 * \param[in] absoluteTolerance Absolute tolerance for comparison.
 */
void compareMtop(FILE* fp, const gmx_mtop_t& mtop1, const gmx_mtop_t& mtop2, real relativeTolerance, real absoluteTolerance);

/*! \brief Check perturbation parameters in topology.
 *
 * \param[in] fp File pointer to write to.
 * \param[in] mtop1 Topology to check perturbation parameters in.
 * \param[in] relativeTolerance Relative tolerance for comparison.
 * \param[in] absoluteTolerance Absolute tolerance for comparison.
 */
void compareMtopAB(FILE* fp, const gmx_mtop_t& mtop1, real relativeTolerance, real absoluteTolerance);

/*! \brief Compare groups.
 *
 * \param[in] fp File pointer to write to.
 * \param[in] g0 First group for comparison.
 * \param[in] g1 Second group for comparison.
 * \param[in] natoms0 Number of atoms for first group.
 * \param[in] natoms1 Number of atoms for second group.
 */
void compareAtomGroups(FILE* fp, const SimulationGroups& g0, const SimulationGroups& g1, int natoms0, int natoms1);

//! Typedef for gmx_localtop in analysis tools.
using ExpandedTopologyPtr = std::unique_ptr<gmx_localtop_t>;

void copy_moltype(const gmx_moltype_t* src, gmx_moltype_t* dst);

#endif
