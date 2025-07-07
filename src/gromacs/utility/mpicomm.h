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
/*! \libinternal \file
 * \brief Declares the MpiComm class
 *
 * This class refers to an MPI communicator and derived data
 * and provides summation reductions that can be hierarchical.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_utility
 * \inlibraryapi
 */
#ifndef GMX_UTILITY_MPICOMM_H
#define GMX_UTILITY_MPICOMM_H

#include <memory>
#include <optional>
#include <string>

#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

template<typename T>
class ArrayRef;

/*! \brief MPI communicator with optimized reduction functionality
 *
 * This class refers to an MPI communicator and derived data
 * and provides summation reductions that can be hierarchical.
 *
 * \note The lifetime of the communicator should be managed externally.
 */
class MpiComm
{
public:
    //! Helper type to differentiate single-rank constructor
    struct SingleRank
    {
    };

    //! Constructs an object for MPI communicator \p comm
    explicit MpiComm(MPI_Comm comm);

    /*! Constructs a communicator with a single rank
     *
     * \note The \p comm() method will return MPI_COMM_NULL. This communicator
     *       should not be passed to MPI calls. Thus any code that might take
     *       an \p MpiComm object constructed by this constructor should check
     *       for \p size() > 1 before calling MPI calls. Note that the
     *       \p sumReduce() methods of this class are safe to call.
     */
    MpiComm(SingleRank /*ignored*/);

    //! Copy constructor
    MpiComm(const MpiComm& other);

    //! Move constructor
    MpiComm(MpiComm&& other) = delete;

    //! Destructor
    ~MpiComm();

    //! Copy assignment
    MpiComm& operator=(const MpiComm& other);

    //! Move assignment
    MpiComm& operator=(MpiComm&& other) noexcept;

    /* \brief ! Returns the MPI communicator
     *
     * \note This return value is MPI_COMM_NULL when \p SingleRank was used to construct
     *       the object. Check for \p size() > 1 before passing such objects to MPI calls.
     */
    MPI_Comm comm() const { return comm_; }

    //! Returns the number of ranks in the communicator
    int size() const { return size_; }

    //! Returns our rank id in the communicator
    int rank() const { return rank_; }

    //! Returns whether our rank is the main rank
    bool isMainRank() const { return rank_ == mainRank_; }

    //! Returns the rank id of the main rank
    int mainRank() const { return mainRank_; }

    /*! \brief Initializes hierarchical reductions
     *
     * When running on system with many nodes, it can be faster to reduce
     * within nodes first and the reduce over nodes using one rank per node.
     * This method initializes such hierarchical reductions. The request is only
     * granted when there is actually such a hierarchy and then is a sufficient
     * number of ranks involved. When that is the case, all later calls to
     * \p sumReduce() will use hierarchical redution. When not, a single reduction
     * will be used.
     *
     * \param physicalNodeIdHash  An integer that should be unique for each physical node
     */
    void initializeHierarchicalReductions(int physicalNodeIdHash);

    //! Returns a string describing the hierarchical reduction, when it is active
    std::optional<std::string> getHierarchicalReductionsReport() const;

    //! Reduce ints over all ranks, results on all ranks, can use hierarchical reduction
    void sumReduce(ArrayRef<int> data) const;

    //! Reduce floats over all ranks, results on all ranks, can use hierarchical reduction
    void sumReduce(ArrayRef<float> data) const;

    //! Reduce doubles over all ranks, results on all ranks, can use hierarchical reduction
    void sumReduce(ArrayRef<double> data) const;

    //! Reduce ints over all ranks, results on all ranks, can use hierarchical reduction
    void sumReduce(std::size_t nr, int* data) const;

    //! Reduce floats over all ranks, results on all ranks, can use hierarchical reduction
    void sumReduce(std::size_t nr, float* data) const;

    //! Reduce doubles over all ranks, results on all ranks, can use hierarchical reduction
    void sumReduce(std::size_t nr, double* data) const;

private:
    //! Machinery for hierarchical reductions over physical nodes
    class HierarchicalReducer;

    //! Rank id of the main rank
    static constexpr int mainRank_ = 0;

    //! Equivalent of MPI_Allreduce in place using summation
    template<typename T>
    void sumReduceTemplated(std::size_t nr, T* data) const;
    //! The MPI communicator, can be MPI_COMM_NULL when there is a single rank
    MPI_Comm comm_ = MPI_COMM_NULL;
    //! The number of ranks in \p comm_
    int size_ = 1;
    //! Our rank id in \p comm_
    int rank_ = 0;
    //! The, optional, machinery for hierarchical reductions
    std::unique_ptr<HierarchicalReducer> hierarchicalReducer_;
};

} // namespace gmx

#endif
