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

#ifndef GMX_NBNXM_PAIRLIST_H
#define GMX_NBNXM_PAIRLIST_H

#include "config.h"

#include <cstddef>

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/defaultinitializationallocator.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "nbnxm_enums.h"

struct t_nblist;

namespace gmx
{

struct NbnxmPairlistCpuWork;
struct NbnxmPairlistGpuWork;

/*! \brief Cache-line protection buffer
 *
 * A buffer data structure of 64 bytes
 * to be placed at the beginning and end of structs
 * to avoid cache invalidation of the real contents
 * of the struct by writes to neighboring memory.
 */
typedef struct
{
    //! Unused field used to create space to protect cache lines that are in use
    int dummy[16];
} gmx_cache_protect_t;

/*! \brief This is the actual cluster-pair list j-entry.
 *
 * cj is the j-cluster.
 * The interaction bits in excl are indexed i-major, j-minor.
 * The cj entries are sorted such that ones with exclusions come first.
 * This means that once a full mask (=NBNXN_INTERACTION_MASK_ALL)
 * is found, all subsequent j-entries in the i-entry also have full masks.
 */
struct nbnxn_cj_t
{
    //! The j-cluster
    int cj;
    //! The exclusion (interaction) bits
    unsigned int excl;
};

//! Simple j-cluster list
class JClusterList
{
public:
    //! The list of packed j-cluster groups
    FastVector<nbnxn_cj_t> list_;

    //! Return the j-cluster index for \c index from the pack list
    int cj(int index) const { return list_[index].cj; }
    //! Return the exclusion mask for \c index
    const unsigned int& excl(int index) const { return list_[index].excl; }
    //! Return the exclusion mask for \c index
    unsigned int& excl(int index) { return list_[index].excl; }
    //! Return the size of the list (not the number of packed elements)
    Index size() const noexcept { return list_.size(); }
    //! Return whether the list is empty
    bool empty() const noexcept { return size() == 0; }
    //! Resize the list
    void resize(Index count) { list_.resize(count); }
    //! Add a new element to the list
    void push_back(const decltype(list_)::value_type& value) { list_.push_back(value); }
};

/*! \brief Constants for interpreting interaction flags
 *
 * In nbnxn_ci_t the integer shift contains the shift in the lower 7 bits.
 * The upper bits contain information for non-bonded kernel optimization.
 * Simply calculating LJ and Coulomb for all pairs in a cluster pair is fine.
 * But three flags can be used to skip interactions, currently only for subc=0
 * !(shift & NBNXN_CI_DO_LJ(subc))   => we can skip LJ for all pairs
 * shift & NBNXN_CI_HALF_LJ(subc)    => we can skip LJ for the second half of i
 * !(shift & NBNXN_CI_DO_COUL(subc)) => we can skip Coulomb for all pairs
 */
//! \{
#define NBNXN_CI_SHIFT 127
#define NBNXN_CI_DO_LJ(subc) (1 << (7 + 3 * (subc)))
#define NBNXN_CI_HALF_LJ(subc) (1 << (8 + 3 * (subc)))
#define NBNXN_CI_DO_COUL(subc) (1 << (9 + 3 * (subc)))
//! \}

/*! \brief Cluster-pair Interaction masks
 *
 * Bit i*j-cluster-size + j tells if atom i and j interact.
 */
//! \{
// TODO: Rename according to convention when moving into Nbnxn namespace
//! All interaction mask is the same for all kernels
constexpr unsigned int NBNXN_INTERACTION_MASK_ALL = 0xffffffffU;
//! \}

/*! \brief Lower limit for square interaction distances in nonbonded kernels.
 *
 * For smaller values we will overflow when calculating r^-1 or r^-12, but
 * to keep it simple we always apply the limit from the tougher r^-12 condition.
 */
#if GMX_DOUBLE
// Some double precision SIMD architectures use single precision in the first
// step, so although the double precision criterion would allow smaller rsq,
// we need to stay in single precision with some margin for the N-R iterations.
constexpr double c_nbnxnMinDistanceSquared = 1.0e-36;
#else
// The worst intermediate value we might evaluate is r^-12, which
// means we should ensure r^2 stays above pow(GMX_FLOAT_MAX,-1.0/6.0)*1.01 (some margin)
constexpr float c_nbnxnMinDistanceSquared = 3.82e-07F; // r > 6.2e-4
#endif

//! Whether we want to use GPU for neighbour list sorting
constexpr bool nbnxmSortListsOnGpu()
{
    return (GMX_GPU_CUDA || GMX_GPU_SYCL);
}

/*! \internal
 * \brief Simple pair-list i-unit
 */
struct nbnxn_ci_t
{
    //! i-cluster
    int ci;
    //! Shift vector index plus possible flags, see above
    int shift;
    //! Start index into cj
    int cj_ind_start;
    //! End index into cj
    int cj_ind_end;
};

//! Grouped pair-list i-unit
struct nbnxn_sci_t
{
    //! Returns the number of j-cluster groups in this entry
    int numJClusterGroups() const { return cjPackedEnd - cjPackedBegin; }

    //! Check if two instances are the same.
    bool operator==(const nbnxn_sci_t& other) const
    {
        return sci == other.sci && shift == other.shift && cjPackedBegin == other.cjPackedBegin
               && cjPackedEnd == other.cjPackedEnd;
    }

    //! i-super-cluster
    int sci;
    //! Shift vector index plus possible flags
    int shift;
    //! Start index into cjPacked
    int cjPackedBegin;
    //! End index into cjPacked (ie. one past the last element)
    int cjPackedEnd;
};

//! Interaction data for a j-group for one warp
struct nbnxn_im_ei_t
{
    //! The i-cluster interactions mask for 1 warp
    unsigned int imask = 0U;
    //! Index into the exclusion array for 1 warp, default index 0 which means no exclusions
    int excl_ind = 0;
    //! Check if two instances are the same.
    bool operator==(const nbnxn_im_ei_t& other) const
    {
        return imask == other.imask && excl_ind == other.excl_ind;
    }
};

//! Packed j-cluster list element
struct nbnxn_cj_packed_t
{
    //! The packed j-clusters
    int cj[c_nbnxnGpuJgroupSize];
    //! The i-cluster mask data for 2 warps
    nbnxn_im_ei_t imei[c_nbnxnGpuClusterpairSplit];
    //! Check if two instances are the same.
    bool operator==(const nbnxn_cj_packed_t& other) const
    {
        return std::equal(std::begin(imei), std::end(imei), std::begin(other.imei), std::end(other.imei))
               && std::equal(std::begin(cj), std::end(cj), std::begin(other.cj), std::end(other.cj));
    }
};

/*! \brief Packed j-cluster list
 *
 * Four j-cluster indices are stored per integer in an nbnxn_cj_packed_t.
 */
class PackedJClusterList
{
public:
    explicit PackedJClusterList(const PinningPolicy pinningPolicy) : list_({}, { pinningPolicy }) {}
    //! The list of packed j-cluster groups
    HostVector<nbnxn_cj_packed_t> list_;
    //! Return the j-cluster index for \c index from the pack list
    int cj(const int index) const
    {
        return list_[index / c_nbnxnGpuJgroupSize].cj[index & (c_nbnxnGpuJgroupSize - 1)];
    }
    //! Return the i-cluster interaction mask for the first cluster in \c index
    unsigned int imask0(const int index) const
    {
        return list_[index / c_nbnxnGpuJgroupSize].imei[0].imask;
    }
    //! Return the size of the list (not the number of packed elements)
    Index size() const noexcept { return list_.size(); }
    //! Return whether the list is empty
    bool empty() const noexcept { return size() == 0; }
    //! Resize the packed list
    void resize(Index count) { list_.resize(count); }
    //! Add a new element to the packed list
    void push_back(const decltype(list_)::value_type& value) { list_.push_back(value); }
};

//! Struct for storing the atom-pair interaction bits for a cluster pair in a GPU pairlist
struct nbnxn_excl_t
{
    //! Constructor, sets no exclusions, so all atom pairs interacting
    MSVC_DIAGNOSTIC_IGNORE(26495) // pair is not being initialized!
    nbnxn_excl_t()
    {
        for (unsigned int& pairEntry : pair)
        {
            pairEntry = NBNXN_INTERACTION_MASK_ALL;
        }
    }
    MSVC_DIAGNOSTIC_RESET

    //! Topology exclusion interaction bits per warp
    unsigned int pair[c_nbnxnGpuExclSize];
    //! Check if two instances are the same.
    bool operator==(const nbnxn_excl_t& other) const
    {
        return std::equal(std::begin(pair), std::end(pair), std::begin(other.pair), std::end(other.pair));
    }
};

//! Cluster pairlist type for use on CPUs
struct NbnxnPairlistCpu
{
    NbnxnPairlistCpu(int iClusterSize);

    //! Cache protection
    gmx_cache_protect_t cp0;

    //! The number of atoms per i-cluster
    int na_ci;
    //! The number of atoms per j-cluster
    int na_cj;
    //! The radius for constructing the list
    real rlist;
    //! The i-cluster list
    FastVector<nbnxn_ci_t> ci;
    //! The outer, unpruned i-cluster list
    FastVector<nbnxn_ci_t> ciOuter;

    //! The j-cluster list
    JClusterList cj;
    //! The outer, unpruned j-cluster list
    FastVector<nbnxn_cj_t> cjOuter;
    //! The number of j-clusters that are used by ci entries in this list, will be <= cj.list.size()
    int ncjInUse;

    //! Working data storage for list construction
    std::unique_ptr<NbnxmPairlistCpuWork> work;

    //! Cache protection
    gmx_cache_protect_t cp1;
};

/* Cluster pairlist type, with extra hierarchies, for on the GPU
 *
 * NOTE: for better performance when combining lists over threads,
 *       all vectors should use default initialization. But when
 *       changing this, excl should be initialized when adding entries.
 */
struct NbnxnPairlistGpu
{
    /*! \brief Constructor
     *
     * \param[in] pinningPolicy  Sets the pinning policy for all buffers used on the GPU
     */
    NbnxnPairlistGpu(PinningPolicy pinningPolicy);

    //! Cache protection
    gmx_cache_protect_t cp0;

    //! The number of atoms per i-cluster
    int na_ci;
    //! The number of atoms per j-cluster
    int na_cj;
    //! The number of atoms per super cluster
    int na_sc;
    //! The radius for constructing the list
    real rlist;
    //! The i-super-cluster list, indexes into cjPacked list;
    HostVector<nbnxn_sci_t> sci;
    //! The list of packed j-cluster groups
    PackedJClusterList cjPacked;
    //! Atom interaction bits (non-exclusions)
    HostVector<nbnxn_excl_t> excl;
    //! The total number of i-clusters
    int nci_tot;

    //! Working data storage for list construction
    std::unique_ptr<NbnxmPairlistGpuWork> work;

    //! Cache protection
    gmx_cache_protect_t cp1;
};

} // namespace gmx

#endif
