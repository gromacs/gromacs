/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
#include "gmxpre.h"

#include "mpicomm.h"

#include "config.h"

#include <limits>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/mpitypes.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

class MpiComm::HierarchicalReducer
{
public:
    /*! \brief
     * Constructor
     *
     * Note that this class takes ownership of the MPI communicators and will free them
     * on destruction.
     *
     * \param[in] commIntra  Communicator covering all ranks within the physical node
     * \param[in] commInter  A communicator between all ranks that have rank 0 in \p commIntra,
     *                       pass MPI_COMM_NULL for ranks that do not have rank 0 in \p commIntra
     * \param[in] numPhysicalNodes  The size of \p commInter
     */
    HierarchicalReducer(MPI_Comm commIntra, MPI_Comm commInter, int numPhysicalNodes);

    //! Constructor, note that this duplicates the communicators
    HierarchicalReducer(const HierarchicalReducer& other);

    //! Destructor, frees the MPI communicators
    ~HierarchicalReducer();

    //! Disallow copy assignment to avoid communicator ownership issues
    HierarchicalReducer& operator=(const HierarchicalReducer&) = delete;

    //! Return a description of the hierarchical setup
    std::string getReport(int numRanksInSim) const;

    //! Equivalent of MPI_Allreduce in place using summation
    template<typename T>
    void sumReduceTemplated(std::size_t nr, T* data) const;

private:
#if GMX_MPI
    //! MPI communicator within a physical node
    MPI_Comm commIntra_;
    //! MPI communicator over nodes, only used between ranks 0 of each physical node
    MPI_Comm commInter_;
#endif
    //! The number of physical nodes
    int numPhysicalNodes_;
};

MpiComm::MpiComm(MPI_Comm comm) : comm_(comm)
{
#if GMX_MPI
    // Handle special case of uninitialized thread-MPI which can occur in tests
    if (!GMX_THREAD_MPI || gmx_mpi_initialized())
    {
        if (comm == MPI_COMM_NULL)
        {
            GMX_THROW(InvalidInputError(
                    "MpiComm(MPI_Comm) can not be called with MPI_COMM_NULL as argument"));
        }

        MPI_Comm_size(comm, &size_);

        if (size_ <= 0)
        {
            GMX_THROW(InvalidInputError("Expect a communicator with at least 1 rank"));
        }

        MPI_Comm_rank(comm, &rank_);
    }
#endif
}

MpiComm::MpiComm(SingleRank /*unused*/) : comm_(MPI_COMM_NULL), size_(1), rank_(0) {}

MpiComm::MpiComm(const MpiComm& other)
{
    *this = other;
}

MpiComm::~MpiComm() = default;

MpiComm& MpiComm::operator=(const MpiComm& other)
{
    comm_ = other.comm_;
    size_ = other.size_;
    rank_ = other.rank_;

    if (&other != this)
    {
        hierarchicalReducer_.reset(nullptr);

        if (other.hierarchicalReducer_)
        {
            hierarchicalReducer_ = std::make_unique<HierarchicalReducer>(*other.hierarchicalReducer_);
        }
    }

    return *this;
}

MpiComm& MpiComm::operator=(MpiComm&& other) noexcept = default;

template<typename T>
void MpiComm::sumReduceTemplated(const std::size_t nr, T* data) const
{
    if (size_ == 1)
    {
        // Only a single rank, no reduction
        return;
    }

#if GMX_MPI
    if (hierarchicalReducer_ == nullptr)
    {
        constexpr std::size_t maxSignedInt = std::numeric_limits<int>::max();

        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Allreduce(MPI_IN_PLACE, data + written, chunk, mpiType<T>(), MPI_SUM, comm_);
            written += chunk;
            remain -= chunk;
        }
    }
    else
    {
        hierarchicalReducer_->sumReduceTemplated(nr, data);
    }
#else
    GMX_RELEASE_ASSERT(false, "No reductions expected without MPI");

    GMX_UNUSED_VALUE(nr);
    GMX_UNUSED_VALUE(data);
#endif
}

template<typename T>
void MpiComm::HierarchicalReducer::sumReduceTemplated(const std::size_t nr, T* data) const
{
#if GMX_MPI
    constexpr std::size_t maxSignedInt = std::numeric_limits<int>::max();

    MPI_Datatype datatype = mpiType<T>();

    for (std::size_t written = 0, remain = nr; remain > 0;)
    {
        const bool isIntraRankZero = (commInter_ != MPI_COMM_NULL);

        std::size_t chunk = std::min(remain, maxSignedInt);
        // Unfortunately MPI requires MPI_IN_PLACE as sendbuf, so we need conditionals ...
        void* sendbuf = (isIntraRankZero ? MPI_IN_PLACE : data + written);
        void* recvbuf = (isIntraRankZero ? data + written : nullptr);
        MPI_Reduce(sendbuf, recvbuf, chunk, datatype, MPI_SUM, 0, commIntra_);
        if (isIntraRankZero)
        {
            /* Sum the roots of the internal (intra) buffers. */
            MPI_Allreduce(MPI_IN_PLACE, data + written, chunk, datatype, MPI_SUM, commInter_);
        }
        written += chunk;
        remain -= chunk;
    }

    for (std::size_t written = 0, remain = nr; remain > 0;)
    {
        std::size_t chunk = std::min(remain, maxSignedInt);
        MPI_Bcast(data + written, chunk, datatype, 0, commIntra_);
        written += chunk;
        remain -= chunk;
    }
#else
    GMX_RELEASE_ASSERT(false, "No reductions expected without MPI");

    GMX_UNUSED_VALUE(nr);
    GMX_UNUSED_VALUE(data);
#endif
}

void MpiComm::sumReduce(ArrayRef<int> data) const
{
    sumReduceTemplated(data.size(), data.data());
}

void MpiComm::sumReduce(ArrayRef<float> data) const
{
    sumReduceTemplated(data.size(), data.data());
}

void MpiComm::sumReduce(ArrayRef<double> data) const
{
    sumReduceTemplated(data.size(), data.data());
}

void MpiComm::sumReduce(std::size_t nr, int* data) const
{
    sumReduceTemplated(nr, data);
}

void MpiComm::sumReduce(std::size_t nr, float* data) const
{
    sumReduceTemplated(nr, data);
}

void MpiComm::sumReduce(std::size_t nr, double* data) const
{
    sumReduceTemplated(nr, data);
}

MpiComm::HierarchicalReducer::HierarchicalReducer(MPI_Comm  commIntra,
                                                  MPI_Comm  commInter,
                                                  const int numPhysicalNodes) :
#if GMX_MPI
    commIntra_(commIntra),
    commInter_(commInter),
#endif
    numPhysicalNodes_(numPhysicalNodes)
{
#if GMX_MPI
    int rank;
    MPI_Comm_rank(commIntra, &rank);
    if ((rank == 0) != (commInter != MPI_COMM_NULL))
    {
        GMX_THROW(
                InvalidInputError("Only rank withs rank 0 in commIntra should have commInter set"));
    }

    if (commInter != MPI_COMM_NULL)
    {
        int numRanks;
        MPI_Comm_size(commInter, &numRanks);
        if (numRanks != numPhysicalNodes)
        {
            GMX_THROW(InvalidInputError("numRanks does not match the size of commInter"));
        }
    }
#else
    GMX_UNUSED_VALUE(commIntra);
    GMX_UNUSED_VALUE(commInter);
#endif
}

MpiComm::HierarchicalReducer::HierarchicalReducer(const HierarchicalReducer& other)
{
#if GMX_MPI
    MPI_Comm_dup(other.commIntra_, &commIntra_);
    if (other.commInter_ != MPI_COMM_NULL)
    {
        MPI_Comm_dup(other.commInter_, &commInter_);
    }
    else
    {
        commInter_ = MPI_COMM_NULL;
    }
#endif
    numPhysicalNodes_ = other.numPhysicalNodes_;
}

MpiComm::HierarchicalReducer::~HierarchicalReducer()
{
#if GMX_MPI
    MPI_Comm_free(&commIntra_);
    if (commInter_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&commInter_);
    }
#endif
}

void MpiComm::initializeHierarchicalReductions(const int physicalNodeIdHash)
{
    if (size() <= 2 || std::getenv("GMX_NO_NODECOMM") != nullptr)
    {
        return;
    }

#if GMX_MPI
    /* Many MPI implementations do not optimize MPI_Allreduce
     * (and probably also other global communication calls)
     * for multi-core nodes connected by a network.
     * We can optimize such communication by using one MPI call
     * within each node and one between the nodes.
     * The most up to date performance check for this hierarchical reduction was
     * performend in 2025 on a HPE Cray EX with Slingshot network using Cray MPICH.
     * On 8 nodes the improvement of the hierarchical reduction is 0.03 ms per call,
     * both with 28 and 14 PP-ranks per node.
     */

    // TODO PhysicalNodeCommunicator could be extended/used to handle
    // the need for per-node per-group communicators.

    /* The intra-node communicator, split on physical node id */
    MPI_Comm commIntra;
    MPI_Comm_split(comm_, physicalNodeIdHash, rank_, &commIntra);
    int rankIntra;
    MPI_Comm_rank(commIntra, &rankIntra);
    /* The inter-node communicator, split on rank_intra.
     * We actually only need the one for rank=0,
     * but it is easier to create them all.
     */
    MPI_Comm commInter;
    MPI_Comm_split(comm_, rankIntra, rank_, &commInter);
    int numRanksInter;
    MPI_Comm_size(commInter, &numRanksInter);

    // Currently only ranks with rankIntra=0 have the correct node count
    int numPhysicalNodes = numRanksInter;
    MPI_Bcast(&numPhysicalNodes, 1, MPI_INT, 0, comm_);

    // We only need the inter-rank communicator on the first rank of each node
    if (rankIntra > 0)
    {
        MPI_Comm_free(&commInter);
    }

    if (numPhysicalNodes > 1 && numPhysicalNodes < size_)
    {
        hierarchicalReducer_ =
                std::make_unique<HierarchicalReducer>(commIntra, commInter, numPhysicalNodes);
    }
    else
    {
        // We don't actually have a hierarchy
        MPI_Comm_free(&commIntra);
        if (commInter != MPI_COMM_NULL)
        {
            MPI_Comm_free(&commInter);
        }
    }
#else
    GMX_UNUSED_VALUE(physicalNodeIdHash);
#endif
}

std::string MpiComm::HierarchicalReducer::getReport(const int numRanksInSim) const
{
    return formatString("Using two-step summing over %d groups of on average %.1f ranks",
                        numPhysicalNodes_,
                        double(numRanksInSim) / double(numPhysicalNodes_));
}

std::optional<std::string> MpiComm::getHierarchicalReductionsReport() const
{
    if (hierarchicalReducer_)
    {
        return hierarchicalReducer_->getReport(size_);
    }
    else
    {
        return {};
    }
}

} // namespace gmx
