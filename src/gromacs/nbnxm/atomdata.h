/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_NBNXN_ATOMDATA_H
#define GMX_NBNXN_ATOMDATA_H

#include <cstdio>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/real.h"

#include "gpu_types.h"

namespace gmx
{
class MDLogger;
}

struct nbnxn_atomdata_t;
struct nonbonded_verlet_t;
struct t_mdatoms;
struct tMPI_Atomic;

class GpuEventSynchronizer;

namespace Nbnxm
{
class GridSet;
enum class KernelType;
} // namespace Nbnxm

/* Convenience type for vector with aligned memory */
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

// Struct that holds force and energy output buffers
struct nbnxn_atomdata_output_t
{
    /* Constructor
     *
     * \param[in] kernelType              Type of non-bonded kernel
     * \param[in] numEnergyGroups         The number of energy groups
     * \param[in] simdEnergyBufferStride  Stride for entries in the energy buffers for SIMD kernels
     * \param[in] pinningPolicy           Sets the pinning policy for all buffers used on the GPU
     */
    nbnxn_atomdata_output_t(Nbnxm::KernelType  kernelType,
                            int                numEnergyGroups,
                            int                simdEnergyBUfferStride,
                            gmx::PinningPolicy pinningPolicy);

    gmx::HostVector<real> f;      // f, size natoms*fstride
    gmx::HostVector<real> fshift; // Shift force array, size SHIFTS*DIM
    gmx::HostVector<real> Vvdw;   // Temporary Van der Waals group energy storage
    gmx::HostVector<real> Vc;     // Temporary Coulomb group energy storage
    AlignedVector<real>   VSvdw;  // Temporary SIMD Van der Waals group energy storage
    AlignedVector<real>   VSc;    // Temporary SIMD Coulomb group energy storage
};

/* Block size in atoms for the non-bonded thread force-buffer reduction,
 * should be a multiple of all cell and x86 SIMD sizes (i.e. 2, 4 and 8).
 * Should be small to reduce the reduction and zeroing cost,
 * but too small will result in overhead.
 * Currently the block size is NBNXN_BUFFERFLAG_SIZE*3*sizeof(real)=192 bytes.
 */
#if GMX_DOUBLE
#    define NBNXN_BUFFERFLAG_SIZE 8
#else
#    define NBNXN_BUFFERFLAG_SIZE 16
#endif

/* We store the reduction flags as gmx_bitmask_t.
 * This limits the number of flags to BITMASK_SIZE.
 */
#define NBNXN_BUFFERFLAG_MAX_THREADS (BITMASK_SIZE)

/* Flags for telling if threads write to force output buffers */
typedef struct
{
    int            nflag;       /* The number of flag blocks                         */
    gmx_bitmask_t* flag;        /* Bit i is set when thread i writes to a cell-block */
    int            flag_nalloc; /* Allocation size of cxy_flag                       */
} nbnxn_buffer_flags_t;

/* LJ combination rules: geometric, Lorentz-Berthelot, none */
enum
{
    ljcrGEOM,
    ljcrLB,
    ljcrNONE,
    ljcrNR
};

/* Struct that stores atom related data for the nbnxn module
 *
 * Note: performance would improve slightly when all std::vector containers
 *       in this struct would not initialize during resize().
 */
struct nbnxn_atomdata_t
{ //NOLINT(clang-analyzer-optin.performance.Padding)
    struct Params
    {
        /* Constructor
         *
         * \param[in] pinningPolicy  Sets the pinning policy for all data that might be transfered to a GPU
         */
        Params(gmx::PinningPolicy pinningPolicy);

        // The number of different atom types
        int numTypes;
        // Lennard-Jone 6*C6 and 12*C12 parameters, size numTypes*2*2
        gmx::HostVector<real> nbfp;
        // Combination rule, see enum defined above
        int comb_rule;
        // LJ parameters per atom type, size numTypes*2
        gmx::HostVector<real> nbfp_comb;
        // As nbfp, but with a stride for the present SIMD architecture
        AlignedVector<real> nbfp_aligned;
        // Atom types per atom
        gmx::HostVector<int> type;
        // LJ parameters per atom for fast SIMD loading
        gmx::HostVector<real> lj_comb;
        // Charges per atom, not set with format nbatXYZQ
        gmx::HostVector<real> q;
        // The number of energy groups
        int nenergrp;
        // 2log(nenergrp)
        int neg_2log;
        // The energy groups, one int entry per cluster, only set when needed
        gmx::HostVector<int> energrp;
    };

    // Diagonal and topology exclusion helper data for all SIMD kernels
    struct SimdMasks
    {
        SimdMasks();

        // Helper data for setting up diagonal exclusion masks in the SIMD 4xN kernels
        AlignedVector<real> diagonal_4xn_j_minus_i;
        // Helper data for setting up diaginal exclusion masks in the SIMD 2xNN kernels
        AlignedVector<real> diagonal_2xnn_j_minus_i;
        // Filters for topology exclusion masks for the SIMD kernels
        AlignedVector<uint32_t> exclusion_filter;
        // Filters for topology exclusion masks for double SIMD kernels without SIMD int32 logical support
        AlignedVector<uint64_t> exclusion_filter64;
        // Array of masks needed for exclusions
        AlignedVector<real> interaction_array;
    };

    /* Constructor
     *
     * \param[in] pinningPolicy  Sets the pinning policy for all data that might be transfered to a GPU
     */
    nbnxn_atomdata_t(gmx::PinningPolicy pinningPolicy);

    /* Returns a const reference to the parameters */
    const Params& params() const { return params_; }

    /* Returns a non-const reference to the parameters */
    Params& paramsDeprecated() { return params_; }

    /* Returns the current total number of atoms stored */
    int numAtoms() const { return numAtoms_; }

    /* Return the coordinate buffer, and q with xFormat==nbatXYZQ */
    gmx::ArrayRef<const real> x() const { return x_; }

    /* Return the coordinate buffer, and q with xFormat==nbatXYZQ */
    gmx::ArrayRef<real> x() { return x_; }

    /* Resizes the coordinate buffer and sets the number of atoms */
    void resizeCoordinateBuffer(int numAtoms);

    /* Resizes the force buffers for the current number of atoms */
    void resizeForceBuffers();

private:
    // The LJ and charge parameters
    Params params_;
    // The total number of atoms currently stored
    int numAtoms_;

public:
    int                        natoms_local; /* Number of local atoms                           */
    int                        XFormat;     /* The format of x (and q), enum                      */
    int                        FFormat;     /* The format of f, enum                              */
    gmx_bool                   bDynamicBox; /* Do we need to update shift_vec every step?    */
    gmx::HostVector<gmx::RVec> shift_vec;   /* Shift vectors, copied from t_forcerec              */
    int                        xstride;     /* stride for a coordinate in x (usually 3 or 4)      */
    int                        fstride;     /* stride for a coordinate in f (usually 3 or 4)      */
private:
    gmx::HostVector<real> x_; /* x and possibly q, size natoms*xstride              */

public:
    // Masks for handling exclusions in the SIMD kernels
    const SimdMasks simdMasks;

    /* Output data */
    std::vector<nbnxn_atomdata_output_t> out; /* Output data structures, 1 per thread */

    /* Reduction related data */
    gmx_bool             bUseBufferFlags; /* Use the flags or operate on all atoms     */
    nbnxn_buffer_flags_t buffer_flags;    /* Flags for buffer zeroing+reduc.  */
    gmx_bool             bUseTreeReduce;  /* Use tree for force reduction */
    tMPI_Atomic*         syncStep;        /* Synchronization step for tree reduce */
};

/* Copy na rvec elements from x to xnb using nbatFormat, start dest a0,
 * and fills up to na_round with coordinates that are far away.
 */
void copy_rvec_to_nbat_real(const int* a, int na, int na_round, const rvec* x, int nbatFormat, real* xnb, int a0);

enum
{
    enbnxninitcombruleDETECT,
    enbnxninitcombruleGEOM,
    enbnxninitcombruleLB,
    enbnxninitcombruleNONE
};

/* Initialize the non-bonded atom data structure.
 * The enum for nbatXFormat is in the file defining nbnxn_atomdata_t.
 * Copy the ntypes*ntypes*2 sized nbfp non-bonded parameter list
 * to the atom data structure.
 * enbnxninitcombrule sets what combination rule data gets stored in nbat.
 */
void nbnxn_atomdata_init(const gmx::MDLogger& mdlog,
                         nbnxn_atomdata_t*    nbat,
                         Nbnxm::KernelType    kernelType,
                         int                  enbnxninitcombrule,
                         int                  ntype,
                         const real*          nbfp,
                         int                  n_energygroups,
                         int                  nout);

void nbnxn_atomdata_set(nbnxn_atomdata_t*     nbat,
                        const Nbnxm::GridSet& gridSet,
                        const t_mdatoms*      mdatoms,
                        const int*            atinfo);

/* Copy the shift vectors to nbat */
void nbnxn_atomdata_copy_shiftvec(gmx_bool dynamic_box, rvec* shift_vec, nbnxn_atomdata_t* nbat);

/*! \brief Transform coordinates to xbat layout
 *
 * Creates a copy of the coordinates buffer using short-range ordering.
 *
 * \param[in] gridSet      The grids data.
 * \param[in] locality     If the transformation should be applied to local or non local coordinates.
 * \param[in] fillLocal    Tells if the local filler particle coordinates should be zeroed.
 * \param[in] coordinates  Coordinates in plain rvec format.
 * \param[in,out] nbat     Data in NBNXM format, used for mapping formats and to locate the output buffer.
 */
void nbnxn_atomdata_copy_x_to_nbat_x(const Nbnxm::GridSet& gridSet,
                                     gmx::AtomLocality     locality,
                                     bool                  fillLocal,
                                     const rvec*           coordinates,
                                     nbnxn_atomdata_t*     nbat);

/*! \brief Transform coordinates to xbat layout on GPU
 *
 * Creates a GPU copy of the coordinates buffer using short-range ordering.
 * As input, uses coordinates in plain rvec format in GPU memory.
 *
 * \param[in]     gridSet    The grids data.
 * \param[in]     locality   If the transformation should be applied to local or non local coordinates.
 * \param[in]     fillLocal  Tells if the local filler particle coordinates should be zeroed.
 * \param[in,out] gpu_nbv    The NBNXM GPU data structure.
 * \param[in]     d_x        Coordinates to be copied (in plain rvec format).
 * \param[in]     xReadyOnDevice   Event synchronizer indicating that the coordinates are ready in the device memory.
 */
void nbnxn_atomdata_x_to_nbat_x_gpu(const Nbnxm::GridSet& gridSet,
                                    gmx::AtomLocality     locality,
                                    bool                  fillLocal,
                                    gmx_nbnxn_gpu_t*      gpu_nbv,
                                    DeviceBuffer<float>   d_x,
                                    GpuEventSynchronizer* xReadyOnDevice);

/*! \brief Add the computed forces to \p f, an internal reduction might be performed as well
 *
 * \param[in]  nbat        Atom data in NBNXM format.
 * \param[in]  locality    If the reduction should be performed on local or non-local atoms.
 * \param[in]  gridSet     The grids data.
 * \param[out] totalForce  Buffer to accumulate resulting force
 */
void reduceForces(nbnxn_atomdata_t* nbat, gmx::AtomLocality locality, const Nbnxm::GridSet& gridSet, rvec* totalForce);

/*! \brief Reduce forces on the GPU
 *
 * \param[in]  locality             If the reduction should be performed on local or non-local atoms.
 * \param[out] totalForcesDevice    Device buffer to accumulate resulting force.
 * \param[in]  gridSet              The grids data.
 * \param[in]  pmeForcesDevice      Device buffer with PME forces.
 * \param[in]  dependencyList       List of synchronizers that represent the dependencies the reduction task needs to sync on.
 * \param[in]  gpu_nbv              The NBNXM GPU data structure.
 * \param[in]  useGpuFPmeReduction  Whether PME forces should be added.
 * \param[in]  accumulateForce      Whether there are usefull data already in the total force buffer.
 */
void reduceForcesGpu(gmx::AtomLocality                          locality,
                     DeviceBuffer<float>                        totalForcesDevice,
                     const Nbnxm::GridSet&                      gridSet,
                     void*                                      pmeForcesDevice,
                     gmx::ArrayRef<GpuEventSynchronizer* const> dependencyList,
                     gmx_nbnxn_gpu_t*                           gpu_nbv,
                     bool                                       useGpuFPmeReduction,
                     bool                                       accumulateForce);

/* Add the fshift force stored in nbat to fshift */
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t& nbat, gmx::ArrayRef<gmx::RVec> fshift);

/* Get the atom start index and number of atoms for a given locality */
void nbnxn_get_atom_range(gmx::AtomLocality     atomLocality,
                          const Nbnxm::GridSet& gridSet,
                          int*                  atomStart,
                          int*                  nAtoms);
#endif
