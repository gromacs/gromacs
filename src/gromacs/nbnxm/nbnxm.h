/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

// FIXME: remove the "__" prefix in front of the group def when we move the
//        nonbonded code into separate dir.

/*! \libinternal \defgroup __module_nbnxm Short-range non-bonded interaction module
 * \ingroup group_mdrun
 *
 * \brief Computes forces and energies for short-range pair-interactions
 * based on the Verlet algorithm. The algorithm uses pair-lists generated
 * at fixed intervals as well as various flavors of pair interaction kernels
 * implemented for a wide range of CPU and GPU architectures.
 *
 * The module includes support for flavors of Coulomb and Lennard-Jones interaction
 * treatment implemented for a large range of SIMD instruction sets for CPU
 * architectures as well as in CUDA and OpenCL for GPU architectures.
 * Additionally there is a reference CPU non-SIMD and a reference CPU
 * for GPU pair-list setup interaction kernel.
 *
 * The implementation of the kernels is based on the cluster non-bonded algorithm
 * which in the code is referred to as the NxM algorithms ("nbnxm_" prefix);
 * for details of the algorithm see DOI:10.1016/j.cpc.2013.06.003.
 *
 * Algorithmically, the non-bonded computation has two different modes:
 * A "classical" mode: generate a list every nstlist steps containing at least
 * all atom pairs up to a distance of rlistOuter and compute pair interactions
 * for all pairs that are within the interaction cut-off.
 * A "dynamic pruning" mode: generate an "outer-list" up to cut-off rlistOuter
 * every nstlist steps and prune the outer-list using a cut-off of rlistInner
 * every nstlistPrune steps to obtain a, smaller, "inner-list". This
 * results in fewer interaction computations and allows for a larger nstlist.
 * On a GPU, this dynamic pruning is performed in a rolling fashion, pruning
 * only a sub-part of the list each (second) step. This way it can often
 * overlap with integration and constraints on the CPU.
 * Currently a simple heuristic determines which mode will be used.
 *
 * TODO: add a summary list and brief descriptions of the different submodules:
 * search, CPU kernels, GPU glue code + kernels.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Anca Hamuraru <anca@streamcomputing.eu>
 * \author Teemu Virolainen <teemu@streamcomputing.eu>
 * \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *
 * TODO: add more authors!
 */

/*! \libinternal
 * \defgroup module_nbnxm Non-bonded pair interactions
 * \ingroup group_mdrun
 * \brief
 * Implements non-bonded pair interaction functionality for NxM atom clusters.
 *
 * This module provides methods to, very efficiently, compute non-bonded
 * pair interactions on CPUs as well as accelerators. It also provides
 * a method to construct the NxM atom-cluster pair-list required for
 * computing these non-bonded iteractions.
 */

/*! \libinternal \file
 *
 * \brief This file contains the public interface of the nbnxm module
 * that implements the NxM atom cluster non-bonded algorithm to efficiently
 * compute pair forces.
 *
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_nbnxm
 */


#ifndef GMX_NBNXM_NBNXM_H
#define GMX_NBNXM_NBNXM_H

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "nbnxm_enums.h"

struct DeviceInformation;
class ExclusionChecker;
struct gmx_enerdata_t;
struct gmx_hw_info_t;
struct gmx_mtop_t;
struct gmx_wallcycle;
struct interaction_const_t;
enum class LJCombinationRule;
struct t_commrec;
struct t_nrnb;
struct t_forcerec;
struct t_inputrec;
struct gmx_grppairener_t;
class GpuEventSynchronizer;

namespace gmx
{
class FreeEnergyDispatch;
struct NbnxmGpu;
struct nbnxn_atomdata_t;
class PairSearch;
class PairlistSets;
template<typename>
class ArrayRefWithPadding;
class DeviceStreamManager;
class DomdecZones;
class ForceWithShiftForces;
class ListedForcesGpu;
template<typename>
class ListOfLists;
class MDLogger;
class ObservablesReducerBuilder;
template<typename>
class Range;
class StepWorkload;
class UpdateGroupsCog;

/* \brief The non-bonded setup, also affects the pairlist construction kernel */
struct NbnxmKernelSetup
{
    //! The non-bonded type, also affects the pairlist construction kernel
    NbnxmKernelType kernelType = NbnxmKernelType::NotSet;
    //! Ewald exclusion computation handling type, currently only used for CPU
    EwaldExclusionType ewaldExclusionType = EwaldExclusionType::NotSet;
};

/*! \brief Flag to tell the nonbonded kernels whether to clear the force output buffers */
enum
{
    enbvClearFNo,
    enbvClearFYes
};

/*! \libinternal
 *  \brief Top-level non-bonded data structure for the Verlet-type cut-off scheme. */
struct nonbonded_verlet_t
{
public:
    /*! \brief Constructs an object from its components
     *
     * \param[in] pairlistSets  The pairlist sets, is consumed
     * \param[in] pairSearch    The pairsearch setup, is consumed
     * \param[in] nbat          The atom data, is consumed
     * \param[in] kernelSetup   The non-bonded kernel setup
     * \param[in] exclusionChecker  The FEP exclusion checker, is consumed, can be nullptr
     * \param[in] gpu_nbv       The GPU non-bonded setup, ownership is transferred, can be nullptr
     * \param[in] wcycle        Pointer to wallcycle counters, can be nullptr
     */
    nonbonded_verlet_t(std::unique_ptr<PairlistSets>     pairlistSets,
                       std::unique_ptr<PairSearch>       pairSearch,
                       std::unique_ptr<nbnxn_atomdata_t> nbat,
                       const NbnxmKernelSetup&           kernelSetup,
                       std::unique_ptr<ExclusionChecker> exclusionChecker,
                       NbnxmGpu*                         gpu_nbv,
                       gmx_wallcycle*                    wcycle);

    /*! \brief Constructs an object from its, minimal, components
     *
     * \param[in] pairlistSets  The pairlist sets, is consumed
     * \param[in] pairSearch    The pairsearch setup, is consumed
     * \param[in] nbat          The atom data, is consumed
     * \param[in] kernelSetup   The non-bonded kernel setup
     * \param[in] gpu_nbv       The GPU non-bonded setup, ownership is transferred, can be nullptr
     */
    nonbonded_verlet_t(std::unique_ptr<PairlistSets>     pairlistSets,
                       std::unique_ptr<PairSearch>       pairSearch,
                       std::unique_ptr<nbnxn_atomdata_t> nbat,
                       const NbnxmKernelSetup&           kernelSetup,
                       NbnxmGpu*                         gpu_nbv);

    ~nonbonded_verlet_t();

    //! Returns whether a GPU is use for the non-bonded calculations
    bool useGpu() const { return isGpuKernelType(kernelSetup_.kernelType); }

    //! Returns whether a GPU is emulated for the non-bonded calculations
    bool emulateGpu() const { return kernelSetup_.kernelType == NbnxmKernelType::Cpu8x8x8_PlainC; }

    //! Return whether the pairlist is of simple, CPU type
    bool pairlistIsSimple() const { return !useGpu() && !emulateGpu(); }

    /*! \brief Put the atoms on the pair search grid.
     *
     * Only atoms with indices wihtin \p atomRange in x are put on the grid.
     * When \p updateGroupsCog != nullptr, atoms are put on the grid
     * based on the center of geometry of the group they belong to.
     * Atoms or COGs of groups should be within the bounding box provided,
     * this is checked in debug builds when not using update groups.
     * The atom density is used to determine the grid size when \p gridIndex = 0.
     * When \p atomDensity <= 0, the density is determined from atomEnd-atomStart
     * and the bounding box corners.
     * With domain decomposition, part of the atoms might have migrated,
     * but have not been removed yet. This count is given by \p numAtomsMoved.
     * When \p move[i] < 0 particle i has migrated and will not be put on the grid.
     *
     * \param[in] box           Box used for periodic distance calculations
     * \param[in] gridIndex     The index of the grid to spread to, always 0 except with test
     * particle insertion
     * \param[in] lowerCorner   Atom groups to be gridded should have coordinates >= this corner
     * \param[in] upperCorner   Atom groups to be gridded should have coordinates <= this corner
     * \param[in] updateGroupsCog  Centers of geometry for update groups,
     *                             pass nullptr when not using update groups
     * \param[in] atomRange     Range of atoms to grid, can include atoms moved to other domains
     * \param[in] numGridAtoms  The number of atoms in \p atomRange excluding moved atoms
     * \param[in] atomDensity   An estimate of the atom density, used for performance optimization,
     *                          only used with \p gridIndex = 0
     * \param[in] atomInfo      Atom information flags
     * \param[in] x             Coordinates for atoms to grid
     * \param[in] move          Move flags for atoms, pass nullptr without DD
     */
    void putAtomsOnGrid(const matrix            box,
                        int                     gridIndex,
                        const RVec&             lowerCorner,
                        const RVec&             upperCorner,
                        const UpdateGroupsCog*  updateGroupsCog,
                        Range<int>              atomRange,
                        int                     numGridAtoms,
                        real                    atomDensity,
                        ArrayRef<const int32_t> atomInfo,
                        ArrayRef<const RVec>    x,
                        const int*              move);

    //! Returns the order of the local atoms on the grid
    ArrayRef<const int> getLocalAtomOrder() const;

    //! Sets the order of the local atoms to the order grid atom ordering
    void setLocalAtomOrder() const;

    //! Returns the index position of the atoms on the search grid
    ArrayRef<const int> getGridIndices() const;

    /*! \brief Constructs the pairlist for the given locality
     *
     * When there are no non-self exclusions, \p exclusions can be empty.
     * Otherwise the number of lists in \p exclusions should match the number
     * of atoms when not using DD, or the total number of atoms in the i-zones
     * when using DD.
     *
     * \param[in] iLocality   The interaction locality: local or non-local
     * \param[in] exclusions  Lists of exclusions for every atom.
     * \param[in] step        Used to set the list creation step
     * \param[in,out] nrnb    Flop accounting struct, can be nullptr
     */
    void constructPairlist(InteractionLocality     iLocality,
                           const ListOfLists<int>& exclusions,
                           int64_t                 step,
                           t_nrnb*                 nrnb) const;

    //! Updates all the atom properties in Nbnxm
    void setAtomProperties(ArrayRef<const int>     atomTypes,
                           ArrayRef<const real>    atomCharges,
                           ArrayRef<const int32_t> atomInfo) const;

    /*!\brief Convert the coordinates to NBNXM format for the given locality.
     *
     * The API function for the transformation of the coordinates from one layout to another.
     *
     * \param[in] locality     Whether coordinates for local or non-local atoms should be
     *                         transformed.
     * \param[in] coordinates  Coordinates in plain rvec format to be transformed.
     */
    void convertCoordinates(AtomLocality locality, ArrayRef<const RVec> coordinates);

    /*!\brief Convert the coordinates to NBNXM format on the GPU for the given locality
     *
     * The API function for the transformation of the coordinates from one layout to another in the GPU memory.
     *
     * \param[in] locality        Whether coordinates for local or non-local atoms should be transformed.
     * \param[in] d_x             GPU coordinates buffer in plain rvec format to be transformed.
     * \param[in] xReadyOnDevice  Event synchronizer indicating that the coordinates are ready in the device memory.
     */
    void convertCoordinatesGpu(AtomLocality locality, DeviceBuffer<RVec> d_x, GpuEventSynchronizer* xReadyOnDevice);

    //! Init for GPU version of setup coordinates in Nbnxm
    void atomdata_init_copy_x_to_nbat_x_gpu() const;

    //! Returns a reference to the pairlist sets
    const PairlistSets& pairlistSets() const { return *pairlistSets_; }

    //! Returns whether step is a dynamic list pruning step, for CPU lists
    bool isDynamicPruningStepCpu(int64_t step) const;

    //! Returns whether step is a dynamic list pruning step, for GPU lists
    bool isDynamicPruningStepGpu(int64_t step) const;

    //! Dispatches the dynamic pruning kernel for the given locality, for CPU lists
    void dispatchPruneKernelCpu(InteractionLocality iLocality, ArrayRef<const RVec> shift_vec) const;

    //! Dispatches the dynamic pruning kernel for GPU lists
    void dispatchPruneKernelGpu(int64_t step);

    //! \brief Executes the non-bonded kernel of the GPU or launches it on the GPU
    void dispatchNonbondedKernel(InteractionLocality        iLocality,
                                 const interaction_const_t& ic,
                                 const StepWorkload&        stepWork,
                                 int                        clearF,
                                 ArrayRef<const RVec>       shiftvec,
                                 ArrayRef<real>             repulsionDispersionSR,
                                 ArrayRef<real>             CoulombSR,
                                 t_nrnb*                    nrnb) const;

    //! Executes the non-bonded free-energy kernels, local + non-local, always runs on the CPU
    void dispatchFreeEnergyKernels(const ArrayRefWithPadding<const RVec>& coords,
                                   ForceWithShiftForces*                  forceWithShiftForces,
                                   bool                                   useSimd,
                                   int                                    ntype,
                                   const interaction_const_t&             ic,
                                   ArrayRef<const RVec>                   shiftvec,
                                   ArrayRef<const real>                   nbfp,
                                   ArrayRef<const real>                   nbfp_grid,
                                   ArrayRef<const real>                   chargeA,
                                   ArrayRef<const real>                   chargeB,
                                   ArrayRef<const int>                    typeA,
                                   ArrayRef<const int>                    typeB,
                                   ArrayRef<const real>                   lambda,
                                   gmx_enerdata_t*                        enerd,
                                   const StepWorkload&                    stepWork,
                                   t_nrnb*                                nrnb);

    /*! \brief Add the forces stored in nbat to f, zeros the forces in nbat
     * \param [in] locality         Local or non-local
     * \param [inout] force         Force to be added to
     */
    void atomdata_add_nbat_f_to_f(AtomLocality locality, ArrayRef<RVec> force);

    /*! \brief Get the number of atoms for a given locality
     *
     * \param [in] locality   Local or non-local
     * \returns               The number of atoms for given locality
     */
    int getNumAtoms(AtomLocality locality) const;

    //! Return the kernel setup
    const NbnxmKernelSetup& kernelSetup() const { return kernelSetup_; }

    //! Returns the outer radius for the pair list
    real pairlistInnerRadius() const;

    //! Returns the outer radius for the pair list
    real pairlistOuterRadius() const;

    //! Changes the pair-list outer and inner radius
    void changePairlistRadii(real rlistOuter, real rlistInner) const;

    //! Set up internal flags that indicate what type of short-range work there is.
    void setupGpuShortRangeWork(const ListedForcesGpu* listedForcesGpu, InteractionLocality iLocality) const;

    void setupFepThreadedForceBuffer(int numAtomsForce);

    //! Returns a reference to the nbnxn_atomdata_t object
    nbnxn_atomdata_t& nbat() { return *nbat_; }

    //! Returns a pointer to the NbnxmGpu object, can return nullptr
    const NbnxmGpu* gpuNbv() const { return gpuNbv_; }

    //! Returns a pointer to the NbnxmGpu object, can return nullptr
    NbnxmGpu* gpuNbv() { return gpuNbv_; }

private:
    //! All data related to the pair lists
    std::unique_ptr<PairlistSets> pairlistSets_;
    //! Working data for constructing the pairlists
    std::unique_ptr<PairSearch> pairSearch_;
    //! Atom data
    std::unique_ptr<nbnxn_atomdata_t> nbat_;

    //! The non-bonded setup, also affects the pairlist construction kernel
    NbnxmKernelSetup kernelSetup_;

    //! \brief The non-bonded free-energy kernel dispatcher
    std::unique_ptr<FreeEnergyDispatch> freeEnergyDispatch_;

    //! \brief Checker for exclusions of perturbed pairs
    std::unique_ptr<ExclusionChecker> exclusionChecker_;

    //! \brief Pointer to wallcycle structure.
    gmx_wallcycle* wcycle_;

    //! GPU Nbnxm data, only used with a physical GPU (TODO: use unique_ptr)
    NbnxmGpu* gpuNbv_;
};

/*! \brief Creates an Nbnxm object */
std::unique_ptr<nonbonded_verlet_t> init_nb_verlet(const MDLogger&            mdlog,
                                                   const t_inputrec&          inputrec,
                                                   const t_forcerec&          forcerec,
                                                   const t_commrec*           commrec,
                                                   const gmx_hw_info_t&       hardwareInfo,
                                                   bool                       useGpuForNonbonded,
                                                   const DeviceStreamManager* deviceStreamManager,
                                                   const gmx_mtop_t&          mtop,
                                                   ObservablesReducerBuilder* observablesReducerBuilder,
                                                   ArrayRef<const RVec>       coordinates,
                                                   matrix                     box,
                                                   gmx_wallcycle*             wcycle);

/*! \brief As nbnxn_put_on_grid, but for the non-local atoms
 *
 * with domain decomposition. Should be called after calling
 * nbnxn_search_put_on_grid for the local atoms / home zone.
 */
void nbnxn_put_on_grid_nonlocal(nonbonded_verlet_t*     nb_verlet,
                                const gmx::DomdecZones& zones,
                                ArrayRef<const int32_t> atomInfo,
                                ArrayRef<const RVec>    x);

/*! \brief Check if GROMACS has been built with GPU support.
 *
 * \param[in] error Pointer to error string or nullptr.
 * \todo Move this to NB module once it exists.
 */
bool buildSupportsNonbondedOnGpu(std::string* error);

} // namespace gmx

#endif // GMX_NBNXM_NBNXM_H
