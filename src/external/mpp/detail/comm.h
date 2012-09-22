/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
#pragma once

#include <memory>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <boost/scoped_array.hpp>
#include "message.h"
#include "decls.h"

namespace mpi {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// comm: is the abstraction of the MPI_Comm class
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class comm{
    MPI_Comm m_comm;
    bool     m_initialized;
    int      m_comm_size;
    int      m_rank;

    // Check whether MPI_Init has been called
    inline void check_init() {

        if (m_initialized) { return; }

        int flag;
        MPI_Initialized(&flag);
        assert(flag != 0 &&
            "FATAL: MPI environment not initialized (MPI_Init not called)");

        m_initialized = true;
        MPI_Comm_size(m_comm, &m_comm_size);
        MPI_Comm_rank(m_comm, &m_rank);
    }

    // Check if any of the selections don't make sense
    inline bool check_selection (const std::vector<int>& sel, const int vecSize) {
        for (unsigned int i=0; i<sel.size(); i++) {
            int s = sel[i];
            if (s > vecSize || s < 0) {
                std::cerr << "ERROR in MPI rank '" << rank()
                          << "': Attempted to select an item outside range of vector!\n"
                          << "Vector has a size of '" << vecSize
                          << "' but tried to acces index '" << s << "'\n";
                return false;
            }
        }
        return true;
    }
public:
    comm(MPI_Comm comm):
        m_comm(comm),
        m_initialized(false),
        m_comm_size(-1),
        m_rank(-1) { }

    // MPI_COMM_WORLD
    static comm world;

    inline int rank() {
        check_init();
        return m_rank;
    }

    inline int rank() const {
        assert(m_initialized && "MPI communicator not initialized");
        return m_rank;
    }

    inline int size() {
        check_init();
        return m_comm_size;
    }

    inline int size() const {
        assert(m_initialized && "MPI communicator not initialized");
        return m_comm_size;
    }

    inline bool isMaster() {
        check_init();

        if (rank()==0) { return true; }
        else           { return false; }
    }

    //Converts comm class into an MPI_Comm
    operator MPI_Comm() {return m_comm;}

    inline endpoint operator()( const int& rank_id );

    //Later: Still need to make both send and receive functions
    template <class RawType>
    inline int send (RawType& buffer, int destination, int tag=0);

    template <class RawType>
    inline int recv (RawType& buffer, int source, int tag=0);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Gather functionality for MPP types (types with mpi_type_traits)
//
// SendType should be some MPP type T or an MPP container of MPP type T.
// RecvType should be an MPP container type of the same MPP type T.
// When gathering containers, it is assumed that
//     the receive-buffer-size / comm-size = number of elements to collect
//     from each individual core.
// sendBuffer should be an MPP type that is to be sent.
// recvBuffer should be a container of the same MPP type as what is being sent.
// root is the rank of the core that should receive the message.
// (Containers can be different types, but contained type must be the same.)
// (recvBuffer can be empty on sending-only cores)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template <class SendType, class RecvType>
    inline int gather (SendType& sendBuffer, RecvType& recvBuffer, int root) {
        msg_impl<SendType> m(sendBuffer);
        msg_impl<RecvType> r(recvBuffer);
        datatype<SendType> m_dt(sendBuffer);
        datatype<RecvType> r_dt;
        // gatherv must be used for dynamic types, since each instance is really a different type.
        if (!m.type_static()) {
            std::vector<int> recvCount(size(), recvBuffer.size()/size());
            return gatherv(sendBuffer, recvBuffer, recvCount, root);
        }

#ifndef NDEBUG
        // Check that receive buffer is the correct size.
        if (rank() == root) {
            assert(r.size() == m.size() * size());
            r_dt.reset(recvBuffer);
        }
#endif
        if (MPI_Gather (const_cast<void*>(m.addr()), m.size(), m_dt.get(),
                        const_cast<void*>(r.addr()), m.size(), r_dt.get(), root, *this)
                        != MPI_SUCCESS) {
            std::ostringstream ss;
            ss << "ERROR in MPI rank '" << comm::world.rank()
               << "': Failed to gather message to destination rank '"
               << root << "'";
            throw comm_error( ss.str() );
        }
        return MPI_SUCCESS;
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Gatherv functionality for MPP types (types with mpi_type_traits)
//
// SendType and RecvType should be MPP container types of some MPP type T.
// sendBuffer should be a container of an MPP type where the entire container
//     is to be sent.
// recvBuffer should be a container of the same MPP type as what is being sent.
// recvCount should be a vector with the Nth element representing how many
//     elements the root should expect to receive from the Nth core.
// root is the rank of the core that should receive the message.
// (Containers can be different types, but contained type must be the same.)
// (recvBuffer and recvCount can be empty on sending-only cores)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template <class SendType, class RecvType>
    inline int gatherv (SendType& sendBuffer, RecvType& recvBuffer,
                        const std::vector<int>& recvCount, int root) {
        msg_impl<SendType> m(sendBuffer);
        msg_impl<RecvType> r(recvBuffer);
        size_t c_size = this->size();
        int mpi_status;
        std::vector<int>          m_count (c_size);
        std::vector<int>          r_disps (c_size);
        m_count[root] = m.size();
        datatype<SendType> m_type_root(sendBuffer);
        //TODO: The if statement is too restrictive, a vector<T> where T is static can also be send with this method
        //But at the same time it is also not restrictive enough. Even if m.type_static() is true, this method can't
        //receive into a list.
        if (m.type_static()) {
#ifndef NDEBUG
            // Check that receive buffer is the correct size.
            if (rank() == root) {
                int totalRecvCount = 0;
                for (size_t i=0; i<c_size; i++) {
                    totalRecvCount += recvCount[i];
                }
                assert(r.size() == totalRecvCount * size());
            }
#endif
            //FIXME: this must be wrong. The type needs to be for a single element in the receiving vector
            // not the whole vector.
            datatype<RecvType> r_type_root(recvBuffer);
            for (size_t i=0; i<c_size; i++) {
                r_disps[i+1] = r_disps[i] + recvCount[i];
            }
            mpi_status = MPI_Gatherv (const_cast<void*>(m.addr()),  m_count[root],      m_type_root.get(),
                                      const_cast<void*>(r.addr()), const_cast<int*>(&recvCount.front()),
                                      &r_disps.front(), r_type_root.get(), root, *this);
        } else {
            boost::scoped_array<MPI_Datatype> m_types(new MPI_Datatype[c_size]);
            std::fill(m_types.get(), m_types.get()+c_size, MPI_CHAR); //Using 0 chars for empty messages
            m_types[root] = m_type_root.get();
            boost::scoped_array<datatype_range> r_types(new datatype_range[c_size]);
            std::vector<int> m_disps (c_size);
            std::vector<int> r_count (c_size);
            typedef typename RecvType::iterator RecvIterator;
            RecvIterator typeDisp = recvBuffer.begin();
            if (rank()==root) {
                for (size_t i=0; i<c_size; i++) {
                    r_count[i] = recvCount[i] == 0 ? 0 : 1;
                    RecvIterator typeEnd = typeDisp;
                    std::advance(typeEnd, recvCount[i]);
                    r_types[i].reset(typeDisp, typeEnd, recvBuffer.begin());
                    typeDisp = typeEnd;
                }
            }
            // m_disps and r_disps have all zero values. r_count has all zero values for non-root and for root has all values
            // one or zero (if receive count is zero). m_count is zero for non-root and is the size of the send buffer for root.
            mpi_status = MPI_Alltoallw(const_cast<void*>(m.addr()), &m_count.front(), &m_disps.front(), m_types.get(),
                                       const_cast<void*>(r.addr()), &r_count.front(), &r_disps.front(), const_cast<MPI_Datatype*>(&r_types.get()->get()), *this);
        }
        if (mpi_status == MPI_SUCCESS) {return MPI_SUCCESS;}
        else {
            std::ostringstream ss;
            ss << "ERROR in MPI rank '" << comm::world.rank()
                << "': Failed to gather message to destination rank '"
                << root << "'";
            throw comm_error( ss.str() );
        }
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Gatherv functionality for MPP types (types with mpi_type_traits)
//
// SendType and RecvType should be MPP container types of some MPP type T.
// sendBuffer should be a container of an MPP type that will send either part
//     or the entire container.
// sendSelection refers to which elements in sendBuffer should be sent.
// recvBuffer should be a container of the same MPP type as what is being sent.
// recvCount should be a vector with the Nth element representing how many
//     elements the root should expect to receive from the Nth core.
// recvIndex should be the same size as recvBuffer and should specify where
//     in recvBuffer each element should be placed. It expects that the Nth
//     recvCount specifies how many indices come from the Nth core.
// root is the rank of the core that should receive the message.
// (Containers can be different types, but contained type must be the same.)
// (recvBuffer, recvCount, and recvIndex can be empty on sending-only cores)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //TODO: check whether the two gatherv calls can be combined by using iterators
    template <class SendType, class RecvType>
    inline int gatherv (SendType& sendBuffer, const std::vector<int>& sendIndex,
                        RecvType& recvBuffer, const std::vector<int>& recvCount,
                        const std::vector<int>& recvIndex, int root) {
        msg_impl<SendType> m(sendBuffer);
        msg_impl<RecvType> r(recvBuffer);
        size_t c_size = this->size();
        int mpi_status;
        std::vector<int>          m_count (c_size);
        std::vector<int>          r_count (c_size);
        std::vector<int>          m_disps (c_size);
        std::vector<int>          r_disps (c_size);

        boost::scoped_array<datatype_range> r_types(new datatype_range[c_size]);
        boost::scoped_array<MPI_Datatype> m_types(new MPI_Datatype[c_size]);
        std::fill(m_types.get(), m_types.get()+c_size, MPI_CHAR);

        assert (check_selection (sendIndex, sendBuffer.size()));
        assert (check_selection (recvIndex, recvBuffer.size()));

        datatype_range m_type_root(
                internal::make_permutation_iterator(sendBuffer.begin(),sendIndex.begin()),
                internal::make_permutation_iterator(sendBuffer.begin(),sendIndex.end()),
                sendBuffer.begin());
        m_count[root] = sendIndex.size() == 0 ? 0 : 1;
        m_types[root] = m_type_root.get();

        if (rank()==root) {
            typedef std::vector<int>::const_iterator IdxIterator;
            IdxIterator idxDisp = recvIndex.begin();
            for (size_t i=0; i<c_size; i++) {
                IdxIterator idxEnd = idxDisp;
                std::advance(idxEnd, recvCount[i]); //idxEnd=idxDisp+recvCount[i]
                r_types[i].reset(
                    internal::make_permutation_iterator(recvBuffer.begin(),idxDisp),
                    internal::make_permutation_iterator(recvBuffer.begin(),idxEnd),
                    recvBuffer.begin());
                idxDisp = idxEnd;
                r_count[i] = recvCount[i] == 0 ? 0 : 1;
            }
        }
        //MPI_Gatherv expects that each section of the receive buffer is a multiple of some datatype.
        //(A section is the part coming from one of the ranks). If the container is either not a vector/array
        //or the contained type is not static, it is not possible to construct the sections from one type.
        //And Alltoallw is the only method accepting different datatypes.

        // m_disps and r_disps have all zero values. r_count has all zero values for non-root and for root has all values
        // one or zero (if receive count is zero). m_count is zero for non-root and is the number of elements sent for root.
        mpi_status = MPI_Alltoallw(const_cast<void*>(m.addr()), &m_count.front(), &m_disps.front(), m_types.get(),
                                   const_cast<void*>(r.addr()), &r_count.front(), &r_disps.front(), const_cast<MPI_Datatype*>(&r_types.get()->get()), *this);

        if (mpi_status == MPI_SUCCESS) {return MPI_SUCCESS;}
        else {
            std::ostringstream ss;
            ss << "ERROR in MPI rank '" << comm::world.rank()
               << "': Failed to gather message to destination rank '"
               << root << "'";
            throw comm_error( ss.str() );
        }
    }
};

inline void init(int* argc = 0, char** argv[] = NULL){
    MPI_Init(argc, argv);
    comm::world = comm(MPI_COMM_WORLD);
}
inline void finalize() {
    for (std::vector<MPI_Datatype>::iterator it=cached_types.begin(); it!=cached_types.end(); ++it) {
        MPI_Type_free(&*it);
    }
    MPI_Finalize();
}
const int any = MPI_ANY_SOURCE;

} // end mpi namespace
