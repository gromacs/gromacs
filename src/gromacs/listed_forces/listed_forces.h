/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
/*! \defgroup module_listed_forces Interactions between lists of particles
 * \ingroup group_mdrun
 *
 * \brief Handles computing energies and forces for listed
 * interactions.
 *
 * Located here is the code for
 * - computing energies and forces for interactions between a small
     number of particles, e.g bonds, position restraints and listed
     non-bonded interactions (e.g. 1-4).
 * - high-level functions used by mdrun for computing a set of such
     quantities
 * - managing thread-wise decomposition, thread-local buffer output,
     and reduction of output data across threads.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 *
 */
/*! \libinternal \file
 *
 * \brief This file contains declarations of high-level functions used
 * by mdrun to compute energies and forces for listed interactions.
 *
 * Clients of libgromacs that want to evaluate listed interactions
 * should call functions declared here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 *
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_H
#define GMX_LISTED_FORCES_LISTED_FORCES_H

#include <bitset>
#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/classhelpers.h"

struct bonded_threading_t;
struct gmx_enerdata_t;
struct gmx_ffparams_t;
struct gmx_grppairener_t;
struct gmx_multisim_t;
class history_t;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_nrnb;

namespace gmx
{
class ForceOutputs;
class StepWorkload;
template<typename>
class ArrayRef;
template<typename>
class ArrayRefWithPadding;
} // namespace gmx

/*! \libinternal
 * \brief Class for calculating listed interactions, uses OpenMP parallelization
 *
 * Listed interactions can be divided over multiple instances of ListedForces
 * using the selection flags passed to the constructor.
 */
class ListedForces
{
public:
    //! Enum for selecting groups of listed interaction types
    enum class InteractionGroup : int
    {
        Pairs,     //!< Pair interactions
        Dihedrals, //!< Dihedrals, including cmap
        Angles,    //!< Angles
        Rest,      //!< All listed interactions that are not any of the above
        Count      //!< The number of items above
    };

    //! Type for specifying selections of groups of interaction types
    using InteractionSelection = std::bitset<static_cast<int>(InteractionGroup::Count)>;

    //! Returns a selection with all listed interaction types selected
    static InteractionSelection interactionSelectionAll()
    {
        InteractionSelection is;
        return is.flip();
    }

    /*! \brief Constructor
     *
     * \param[in] ffparams         The force field parameters
     * \param[in] numEnergyGroups  The number of energy groups, used for storage of pair energies
     * \param[in] numThreads       The number of threads used for computed listed interactions
     * \param[in] interactionSelection  Select of interaction groups through bits set
     * \param[in] fplog            Log file for printing env.var. override, can be nullptr
     */
    ListedForces(const gmx_ffparams_t& ffparams,
                 int                   numEnergyGroups,
                 int                   numThreads,
                 InteractionSelection  interactionSelection,
                 FILE*                 fplog);

    //! Move constructor, default, but in the source file to hide implementation classes
    ListedForces(ListedForces&& o) noexcept;

    //! Destructor which is actually default but in the source file to hide implementation classes
    ~ListedForces();

    /*! \brief Copy the listed interactions from \p idef and set up the thread parallelization
     *
     * \param[in] domainIdef     Interaction definitions for all listed interactions to be computed on this domain/rank
     * \param[in] numAtomsForce  Force are, potentially, computed for atoms 0 to \p numAtomsForce
     * \param[in] useGpu         Whether a GPU is used to compute (part of) the listed interactions
     */
    void setup(const InteractionDefinitions& domainIdef, int numAtomsForce, bool useGpu);

    /*! \brief Do all aspects of energy and force calculations for mdrun
     * on the set of listed interactions
     *
     * xWholeMolecules only needs to contain whole molecules when orientation
     * restraints need to be computed and can be empty otherwise.
     */
    void calculate(struct gmx_wallcycle*                     wcycle,
                   const matrix                              box,
                   const t_commrec*                          cr,
                   const gmx_multisim_t*                     ms,
                   gmx::ArrayRefWithPadding<const gmx::RVec> coordinates,
                   gmx::ArrayRef<const gmx::RVec>            xWholeMolecules,
                   t_fcdata*                                 fcdata,
                   const history_t*                          hist,
                   gmx::ForceOutputs*                        forceOutputs,
                   const t_forcerec*                         fr,
                   const struct t_pbc*                       pbc,
                   gmx_enerdata_t*                           enerd,
                   t_nrnb*                                   nrnb,
                   gmx::ArrayRef<const real>                 lambda,
                   gmx::ArrayRef<const real>                 chargeA,
                   gmx::ArrayRef<const real>                 chargeB,
                   gmx::ArrayRef<const bool>                 atomIsPerturbed,
                   gmx::ArrayRef<const unsigned short>       cENER,
                   int                                       nPerturbed,
                   int*                                      global_atom_index,
                   const gmx::StepWorkload&                  stepWork);

    //! Returns whether bonded interactions are assigned to the CPU
    bool haveCpuBondeds() const;

    /*! \brief Returns whether listed forces are computed on the CPU
     *
     * NOTE: the current implementation returns true if there are position restraints
     * or any bonded interactions computed on the CPU.
     */
    bool haveCpuListedForces(const t_fcdata& fcdata) const;

    //! Returns true if there are position, distance or orientation restraints
    bool haveRestraints(const t_fcdata& fcdata) const;

private:
    //! Pointer to the interaction definitions
    InteractionDefinitions const* idef_ = nullptr;
    //! The number of energy groups
    int numEnergyGroups_;
    //! Interaction definitions used for storing selections
    InteractionDefinitions idefSelection_;
    //! Thread parallelization setup, unique_ptr to avoid declaring bonded_threading_t
    std::unique_ptr<bonded_threading_t> threading_;
    //! Tells which interactions to select for computation
    const InteractionSelection interactionSelection_;
    //! Force buffer for free-energy forces
    std::vector<real> forceBufferLambda_;
    //! Shift force buffer for free-energy forces
    std::vector<gmx::RVec> shiftForceBufferLambda_;
    //! Temporary array for storing foreign lambda group pair energies
    std::unique_ptr<gmx_grppairener_t> foreignEnergyGroups_;

    GMX_DISALLOW_COPY_AND_ASSIGN(ListedForces);
};

#endif
