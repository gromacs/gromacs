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
/*! \libinternal \file
 *  \brief
 *  Functionality for per-atom data in the nbnxm module
 *
 *  \author Berk Hess <hess@kth.se>
 *  \ingroup module_nbnxm
 *  \inlibraryapi
 */


#ifndef GMX_NBNXN_ATOMDATA_H
#define GMX_NBNXN_ATOMDATA_H

#include <cstdint>
#include <cstdio>

#include <memory>
#include <optional>
#include <vector>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/real.h"

class GpuEventSynchronizer;

namespace gmx
{

template<bool, bool>
class EnergyAccumulator;
class EnergyGroupsPerCluster;
class MDLogger;
struct NbnxmGpu;
struct nbnxn_atomdata_t;
struct nonbonded_verlet_t;
enum class NbnxmKernelType;
class GridSet;

enum
{
    nbatXYZ,
    nbatXYZQ,
    nbatX4,
    nbatX8
};

//! Stride for coordinate/force arrays with xyz coordinate storage
static constexpr int STRIDE_XYZ = 3;
//! Stride for coordinate/force arrays with xyzq coordinate storage
static constexpr int STRIDE_XYZQ = 4;
//! Size of packs of x, y or z with SIMD 4-grouped packed coordinates/forces
static constexpr int c_packX4 = 4;
//! Size of packs of x, y or z with SIMD 8-grouped packed coordinates/forces
static constexpr int c_packX8 = 8;
//! Stridefor a pack of 4 coordinates/forces
static constexpr int STRIDE_P4 = DIM * c_packX4;
//! Stridefor a pack of 8 coordinates/forces
static constexpr int STRIDE_P8 = DIM * c_packX8;

//! Returns the index in a coordinate array corresponding to atom a
template<int packSize>
static inline int atom_to_x_index(int a)
{
    return DIM * (a & ~(packSize - 1)) + (a & (packSize - 1));
}

/*! \internal
 * \brief Struct that holds force and energy output buffers */
struct nbnxn_atomdata_output_t
{
    /*! \brief Constructor
     *
     * \param[in] kernelType       Type of non-bonded kernel
     * \param[in] numEnergyGroups  The number of energy groups
     * \param[in] pinningPolicy    Sets the pinning policy for all buffers used on the GPU
     */
    nbnxn_atomdata_output_t(NbnxmKernelType kernelType, int numEnergyGroups, PinningPolicy pinningPolicy);

    //! Move constructor
    nbnxn_atomdata_output_t(nbnxn_atomdata_output_t&&) noexcept;

    //! Destructor
    ~nbnxn_atomdata_output_t();


    //! f, size natoms*fstride
    HostVector<real> f;
    //! Shift force array, size c_numShiftVectors*DIM
    HostVector<real> fshift;
    //! Temporary Van der Waals group energy storage
    HostVector<real> Vvdw;
    //! Temporary Coulomb group energy storage
    HostVector<real> Vc;

    //! Accumulator for energy output with a single energy group
    std::unique_ptr<EnergyAccumulator<false, true>> accumulatorSingleEnergies;
    //! Accumulator for energy output with multiple energy groups
    std::unique_ptr<EnergyAccumulator<true, true>> accumulatorGroupEnergies;
};

/*! \brief Block size in atoms for the non-bonded thread force-buffer reduction.
 *
 * Should be a multiple of all cell and x86 SIMD sizes (i.e. 2, 4 and 8).
 * Should be small to reduce the reduction and zeroing cost,
 * but too small will result in overhead.
 * Currently the block size is NBNXN_BUFFERFLAG_SIZE*3*sizeof(real)=192 bytes.
 */
#if GMX_DOUBLE
#    define NBNXN_BUFFERFLAG_SIZE 8
#else
#    define NBNXN_BUFFERFLAG_SIZE 16
#endif

/*! \brief We store the reduction flags as gmx_bitmask_t.
 * This limits the number of flags to BITMASK_SIZE.
 */
#define NBNXN_BUFFERFLAG_MAX_THREADS (BITMASK_SIZE)


//! LJ combination rules
enum class LJCombinationRule : int
{
    //! Geometric
    Geometric,
    //! Lorentz-Berthelot
    LorentzBerthelot,
    //! No rule
    None,
    //! Size of the enum
    Count
};

//! String corresponding to LJ combination rule
const char* enumValueToString(LJCombinationRule enumValue);

/*! \internal
 * \brief Struct that stores atom related data for the nbnxn module
 *
 * Note: performance would improve slightly when all std::vector containers
 *       in this struct would not initialize during resize().
 */
struct nbnxn_atomdata_t
{ //NOLINT(clang-analyzer-optin.performance.Padding)

    nbnxn_atomdata_t(PinningPolicy pinningPolicy);

    ~nbnxn_atomdata_t();

    /*! \internal
     * \brief The actual atom data parameter values */
    struct Params
    {
        /*! \brief Constructor
         *
         * \param[in] pinningPolicy  Sets the pinning policy for all data that might be transfered to a GPU
         */
        Params(PinningPolicy pinningPolicy);

        //! The number of different atom types
        int numTypes;
        //! Lennard-Jone 6*C6 and 12*C12 parameters, size numTypes*2*2
        HostVector<real> nbfp;
        //! Combination rule, see enum defined above
        LJCombinationRule ljCombinationRule;
        //! LJ parameters per atom type, size numTypes*2
        HostVector<real> nbfp_comb;
        //! As nbfp, but with a stride for the present SIMD architecture
        AlignedVector<real> nbfp_aligned;
        //! Atom types per atom
        HostVector<int> type;
        //! LJ parameters per atom for fast SIMD loading
        HostVector<real> lj_comb;
        //! Charges per atom, not set with format nbatXYZQ
        HostVector<real> q;
        //! The number of energy groups
        int numEnergyGroups;
        //! The list of energy groups per i-cluster
        std::unique_ptr<EnergyGroupsPerCluster> energyGroupsPerCluster;
    };

    /*! \internal
     * \brief Diagonal and topology exclusion helper data for all SIMD kernels. */
    struct SimdMasks
    {
        SimdMasks(NbnxmKernelType kernelType);

        //! Helper data for setting up diagonal exclusion masks in the SIMD 4xN kernels
        AlignedVector<real> diagonal_4xn_j_minus_i;
        //! Helper data for setting up diaginal exclusion masks in the SIMD 2xNN kernels
        AlignedVector<real> diagonal_2xnn_j_minus_i;
        //! Filters for topology exclusion masks for the SIMD kernels
        AlignedVector<uint32_t> exclusion_filter;
        //! Filters for topology exclusion masks for double SIMD kernels without SIMD int32 logical support
        AlignedVector<uint64_t> exclusion_filter64;
    };

    /*! \brief Constructor
     *
     * This class only stores one, or no, LJ combination rule parameter list.
     * With LJ-PME the rule for the LJ PME-grid must match the PME grid combination rule
     * and there can be no combination rule for the LJ pair parameters.
     * Without LJ-PME the combination rule should match the combination rule
     * for the LJ parameters.
     * An (release) assertion failure will occur when these conditions are not met.
     *
     * The non-bonded force parameter matrix \p nbfp is a serialized version
     * of a square parameter matrix with a pair of parameters 6*C6, 12*C12 for every
     * atom type pair.
     *
     * \param[in] pinningPolicy      Sets the pinning policy for all data that might be transferred
     *                               to a GPU
     * \param[in] mdlog              The logger
     * \param[in] kernelType         Nonbonded NxN kernel type
     * \param[in] ljCombinationRule  The LJ combination rule parameters to generate,
                                     empty is detect from the LJ parameters
     * \param[in] pmeLJCombinationRule  The LJ combination rule parameters to generate for the LJ PME-grid part
     * \param[in] nbfp               Non-bonded force parameter matrix
     * \param[in] addFillerAtomType  When true, add a filler atom type, when false, \p nbfp should
     *                               have atom-type with index numTypes-1 with all parameters zero
     *                               so that that row and column have only zero values.
     * \param[in] numEnergyGroups    Number of energy groups
     * \param[in] numOutputBuffers   Number of output data structures
     */
    nbnxn_atomdata_t(PinningPolicy                           pinningPolicy,
                     const MDLogger&                         mdlog,
                     NbnxmKernelType                         kernelType,
                     const std::optional<LJCombinationRule>& ljCombinationRule,
                     LJCombinationRule                       pmeLJCombinationRule,
                     ArrayRef<const real>                    nbfp,
                     bool                                    addFillerAtomType,
                     int                                     numEnergyGroups,
                     int                                     numOutputBuffers);

    //! Returns a const reference to the parameters
    const Params& params() const { return params_; }

    //! Returns a non-const reference to the parameters
    Params& paramsDeprecated() { return params_; }

    //! Returns the current total number of atoms stored
    int numAtoms() const { return numAtoms_; }

    //! Returns the number of local atoms
    int numLocalAtoms() const { return numLocalAtoms_; }

    //! Return the coordinate buffer, and q with xFormat==nbatXYZQ
    ArrayRef<const real> x() const { return x_; }

    //! Return the coordinate buffer, and q with xFormat==nbatXYZQ
    ArrayRef<real> x() { return x_; }

    //! Masks for handling exclusions in the SIMD kernels
    const SimdMasks& simdMasks() const { return simdMasks_; }

    /*! \brief Resizes the coordinate buffer and sets the number of atoms
     *
     * \param numAtoms                 The new number of atoms
     * \param domainDecompositionZone  The domain decomposition zone index
     */
    void resizeCoordinateBuffer(int numAtoms, int domainDecompositionZone);

    //! Resizes the force buffers for the current number of atoms
    void resizeForceBuffers();

    //! Returns the output buffer for the given \p thread
    nbnxn_atomdata_output_t& outputBuffer(int thread) { return outputBuffers_[thread]; }

    //! Returns the output buffer for the given \p thread
    const nbnxn_atomdata_output_t& outputBuffer(int thread) const { return outputBuffers_[thread]; }

    //! Returns the list of output buffers
    ArrayRef<const nbnxn_atomdata_output_t> outputBuffers() const { return outputBuffers_; }

    //! Returns whether buffer flags are used
    bool useBufferFlags() const { return useBufferFlags_; }

    //! Returns the vector of buffer flags
    std::vector<gmx_bitmask_t>& bufferFlags() { return bufferFlags_; }

    /*! \brief Add the computed forces to \p f, an internal reduction might be performed as well
     *
     * \param[in]  locality    If the reduction should be performed on local or non-local atoms.
     * \param[in]  gridSet     The grids data.
     * \param[out] totalForce  Buffer to accumulate resulting force
     */
    void reduceForces(AtomLocality locality, const GridSet& gridSet, rvec* totalForce);

    /*! \brief Clears the force buffer.
     *
     * Either the whole buffer is cleared or only the parts used
     * by thread/task \p outputIndex when useBufferFlags() returns true.
     *
     * \param[in]     outputIndex  The index of the output object to clear
     */
    void clearForceBuffer(int outputIndex);

private:
    //! Reduce the output buffers into the first one
    void reduceForcesOverThreads();

    //! The LJ and charge parameters
    Params params_;
    //! The total number of atoms currently stored
    int numAtoms_;
    //! The number of local atoms
    int numLocalAtoms_;

public:
    //! The format of x (and q), enum
    int XFormat;
    //! The format of f, enum
    int FFormat;
    //! Do we need to update shift_vec every step?
    bool bDynamicBox;
    //! Shift vectors, copied from t_forcerec
    HostVector<RVec> shift_vec;
    //! stride for a coordinate in x (usually 3 or 4)
    int xstride;
    //! stride for a coordinate in f (usually 3 or 4)
    int fstride;

private:
    //! x and possibly q, size natoms*xstride
    HostVector<real> x_;

    //! Masks for handling exclusions in the SIMD kernels
    SimdMasks simdMasks_;

    //! Output data structures, 1 per thread
    std::vector<nbnxn_atomdata_output_t> outputBuffers_;

    //! Reduction related data
    //! \{
    //! Use the flags or operate on all atoms
    bool useBufferFlags_;
    //! Flags for buffer zeroing+reduc.
    std::vector<gmx_bitmask_t> bufferFlags_;
    //! \}
};

/*! \brief Copy na rvec elements from x to xnb using nbatFormat, start dest a0,
 * and fills up to na_round with coordinates that are far away.
 */
void copy_rvec_to_nbat_real(const int* a, int na, int na_round, const rvec* x, int nbatFormat, real* xnb, int a0);

//! Sets the atomdata after pair search
void nbnxn_atomdata_set(nbnxn_atomdata_t*       nbat,
                        const GridSet&          gridSet,
                        ArrayRef<const int>     atomTypes,
                        ArrayRef<const real>    atomCharges,
                        ArrayRef<const int32_t> atomInfo);

//! Copy the shift vectors to nbat
void nbnxn_atomdata_copy_shiftvec(bool dynamic_box, ArrayRef<RVec> shift_vec, nbnxn_atomdata_t* nbat);

/*! \brief Transform coordinates to xbat layout
 *
 * Creates a copy of the coordinates buffer using short-range ordering.
 *
 * \param[in] gridSet      The grids data.
 * \param[in] locality     If the transformation should be applied to local or non local coordinates.
 * \param[in] coordinates  Coordinates in plain rvec format.
 * \param[in,out] nbat     Data in NBNXM format, used for mapping formats and to locate the output buffer.
 */
void nbnxn_atomdata_copy_x_to_nbat_x(const GridSet&    gridSet,
                                     AtomLocality      locality,
                                     const rvec*       coordinates,
                                     nbnxn_atomdata_t* nbat);

/*! \brief Transform coordinates to xbat layout on GPU
 *
 * Creates a GPU copy of the coordinates buffer using short-range ordering.
 * As input, uses coordinates in plain rvec format in GPU memory.
 *
 * \param[in]     gridSet    The grids data.
 * \param[in]     locality   If the transformation should be applied to local or non local coordinates.
 * \param[in,out] gpu_nbv    The NBNXM GPU data structure.
 * \param[in]     d_x        Coordinates to be copied (in plain rvec format).
 * \param[in]     xReadyOnDevice   Event synchronizer indicating that the coordinates are ready in the device memory.
 *                                 If there is no need to wait for any event (e.g., the wait has already been
 *                                 enqueued into the appropriate stream), it can be \c nullptr.
 */
void nbnxn_atomdata_x_to_nbat_x_gpu(const GridSet&        gridSet,
                                    AtomLocality          locality,
                                    NbnxmGpu*             gpu_nbv,
                                    DeviceBuffer<RVec>    d_x,
                                    GpuEventSynchronizer* xReadyOnDevice);

//! Add the fshift force stored in nbat to fshift
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t& nbat, ArrayRef<RVec> fshift);

} // namespace gmx

#endif
