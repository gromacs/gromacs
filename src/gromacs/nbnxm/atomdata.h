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

#include <cstdio>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/real.h"

namespace gmx
{
class MDLogger;
}

struct NbnxmGpu;
struct nbnxn_atomdata_t;
struct nonbonded_verlet_t;

class GpuEventSynchronizer;

namespace Nbnxm
{
class GridSet;
enum class KernelType;
} // namespace Nbnxm

//! Convenience type for vector with aligned memory
template<typename T>
using AlignedVector = std::vector<T, gmx::AlignedAllocator<T>>;

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
     * \param[in] kernelType              Type of non-bonded kernel
     * \param[in] numEnergyGroups         The number of energy groups
     * \param[in] simdEnergyBufferStride  Stride for entries in the energy buffers for SIMD kernels
     * \param[in] pinningPolicy           Sets the pinning policy for all buffers used on the GPU
     */
    nbnxn_atomdata_output_t(Nbnxm::KernelType  kernelType,
                            int                numEnergyGroups,
                            int                simdEnergyBufferStride,
                            gmx::PinningPolicy pinningPolicy);

    //! f, size natoms*fstride
    gmx::HostVector<real> f;
    //! Shift force array, size c_numShiftVectors*DIM
    gmx::HostVector<real> fshift;
    //! Temporary Van der Waals group energy storage
    gmx::HostVector<real> Vvdw;
    //! Temporary Coulomb group energy storage
    gmx::HostVector<real> Vc;
    //! Temporary SIMD Van der Waals group energy storage
    AlignedVector<real> VSvdw;
    //! Temporary SIMD Coulomb group energy storage
    AlignedVector<real> VSc;
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
    /*! \internal
     * \brief The actual atom data parameter values */
    struct Params
    {
        /*! \brief Constructor
         *
         * \param[in] pinningPolicy  Sets the pinning policy for all data that might be transfered to a GPU
         */
        Params(gmx::PinningPolicy pinningPolicy);

        //! The number of different atom types
        int numTypes;
        //! Lennard-Jone 6*C6 and 12*C12 parameters, size numTypes*2*2
        gmx::HostVector<real> nbfp;
        //! Combination rule, see enum defined above
        LJCombinationRule ljCombinationRule;
        //! LJ parameters per atom type, size numTypes*2
        gmx::HostVector<real> nbfp_comb;
        //! As nbfp, but with a stride for the present SIMD architecture
        AlignedVector<real> nbfp_aligned;
        //! Atom types per atom
        gmx::HostVector<int> type;
        //! LJ parameters per atom for fast SIMD loading
        gmx::HostVector<real> lj_comb;
        //! Charges per atom, not set with format nbatXYZQ
        gmx::HostVector<real> q;
        //! The number of energy groups
        int nenergrp;
        //! 2log(nenergrp)
        int neg_2log;
        //! The energy groups, one int entry per cluster, only set when needed
        gmx::HostVector<int> energrp;
    };

    /*! \internal
     * \brief Diagonal and topology exclusion helper data for all SIMD kernels. */
    struct SimdMasks
    {
        SimdMasks();

        //! Helper data for setting up diagonal exclusion masks in the SIMD 4xN kernels
        AlignedVector<real> diagonal_4xn_j_minus_i;
        //! Helper data for setting up diaginal exclusion masks in the SIMD 2xNN kernels
        AlignedVector<real> diagonal_2xnn_j_minus_i;
        //! Filters for topology exclusion masks for the SIMD kernels
        AlignedVector<uint32_t> exclusion_filter;
        //! Filters for topology exclusion masks for double SIMD kernels without SIMD int32 logical support
        AlignedVector<uint64_t> exclusion_filter64;
        //! Array of masks needed for exclusions
        AlignedVector<real> interaction_array;
    };

    /*! \brief Constructor
     *
     * \param[in] pinningPolicy      Sets the pinning policy for all data that might be transferred
     *                               to a GPU
     * \param[in] mdlog              The logger
     * \param[in] kernelType         Nonbonded NxN kernel type
     * \param[in] enbnxninitcombrule LJ combination rule
     * \param[in] ntype              Number of atom types
     * \param[in] nbfp               Non-bonded force parameters
     * \param[in] n_energygroups     Number of energy groups
     * \param[in] nout               Number of output data structures
     */
    nbnxn_atomdata_t(gmx::PinningPolicy        pinningPolicy,
                     const gmx::MDLogger&      mdlog,
                     Nbnxm::KernelType         kernelType,
                     int                       enbnxninitcombrule,
                     int                       ntype,
                     gmx::ArrayRef<const real> nbfp,
                     int                       n_energygroups,
                     int                       nout);

    //! Returns a const reference to the parameters
    const Params& params() const { return params_; }

    //! Returns a non-const reference to the parameters
    Params& paramsDeprecated() { return params_; }

    //! Returns the current total number of atoms stored
    int numAtoms() const { return numAtoms_; }

    //! Return the coordinate buffer, and q with xFormat==nbatXYZQ
    gmx::ArrayRef<const real> x() const { return x_; }

    //! Return the coordinate buffer, and q with xFormat==nbatXYZQ
    gmx::ArrayRef<real> x() { return x_; }

    //! Resizes the coordinate buffer and sets the number of atoms
    void resizeCoordinateBuffer(int numAtoms);

    //! Resizes the force buffers for the current number of atoms
    void resizeForceBuffers();

private:
    //! The LJ and charge parameters
    Params params_;
    //! The total number of atoms currently stored
    int numAtoms_;

public:
    //! Number of local atoms
    int natoms_local;
    //! The format of x (and q), enum
    int XFormat;
    //! The format of f, enum
    int FFormat;
    //! Do we need to update shift_vec every step?
    bool bDynamicBox;
    //! Shift vectors, copied from t_forcerec
    gmx::HostVector<gmx::RVec> shift_vec;
    //! stride for a coordinate in x (usually 3 or 4)
    int xstride;
    //! stride for a coordinate in f (usually 3 or 4)
    int fstride;

private:
    //! x and possibly q, size natoms*xstride
    gmx::HostVector<real> x_;

public:
    //! Masks for handling exclusions in the SIMD kernels
    const SimdMasks simdMasks;

    //! Output data structures, 1 per thread
    std::vector<nbnxn_atomdata_output_t> out;

    //! Reduction related data
    //! \{
    //! Use the flags or operate on all atoms
    bool bUseBufferFlags;
    //! Flags for buffer zeroing+reduc.
    std::vector<gmx_bitmask_t> buffer_flags;
    //! \}
};

/*! \brief Copy na rvec elements from x to xnb using nbatFormat, start dest a0,
 * and fills up to na_round with coordinates that are far away.
 */
void copy_rvec_to_nbat_real(const int* a, int na, int na_round, const rvec* x, int nbatFormat, real* xnb, int a0);

//! Describes the combination rule in use by this force field
enum
{
    enbnxninitcombruleDETECT,
    enbnxninitcombruleGEOM,
    enbnxninitcombruleLB,
    enbnxninitcombruleNONE
};

//! Sets the atomdata after pair search
void nbnxn_atomdata_set(nbnxn_atomdata_t*            nbat,
                        const Nbnxm::GridSet&        gridSet,
                        gmx::ArrayRef<const int>     atomTypes,
                        gmx::ArrayRef<const real>    atomCharges,
                        gmx::ArrayRef<const int64_t> atomInfo);

//! Copy the shift vectors to nbat
void nbnxn_atomdata_copy_shiftvec(bool dynamic_box, gmx::ArrayRef<gmx::RVec> shift_vec, nbnxn_atomdata_t* nbat);

/*! \brief Transform coordinates to xbat layout
 *
 * Creates a copy of the coordinates buffer using short-range ordering.
 *
 * \param[in] gridSet      The grids data.
 * \param[in] locality     If the transformation should be applied to local or non local coordinates.
 * \param[in] coordinates  Coordinates in plain rvec format.
 * \param[in,out] nbat     Data in NBNXM format, used for mapping formats and to locate the output buffer.
 */
void nbnxn_atomdata_copy_x_to_nbat_x(const Nbnxm::GridSet& gridSet,
                                     gmx::AtomLocality     locality,
                                     const rvec*           coordinates,
                                     nbnxn_atomdata_t*     nbat);

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
void nbnxn_atomdata_x_to_nbat_x_gpu(const Nbnxm::GridSet&   gridSet,
                                    gmx::AtomLocality       locality,
                                    NbnxmGpu*               gpu_nbv,
                                    DeviceBuffer<gmx::RVec> d_x,
                                    GpuEventSynchronizer*   xReadyOnDevice);

/*! \brief Add the computed forces to \p f, an internal reduction might be performed as well
 *
 * \param[in]  nbat        Atom data in NBNXM format.
 * \param[in]  locality    If the reduction should be performed on local or non-local atoms.
 * \param[in]  gridSet     The grids data.
 * \param[out] totalForce  Buffer to accumulate resulting force
 */
void reduceForces(nbnxn_atomdata_t* nbat, gmx::AtomLocality locality, const Nbnxm::GridSet& gridSet, rvec* totalForce);

//! Add the fshift force stored in nbat to fshift
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t& nbat, gmx::ArrayRef<gmx::RVec> fshift);

//! Get the atom start index and number of atoms for a given locality
void nbnxn_get_atom_range(gmx::AtomLocality     atomLocality,
                          const Nbnxm::GridSet& gridSet,
                          int*                  atomStart,
                          int*                  nAtoms);
#endif
