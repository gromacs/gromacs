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
#include "msg_type.h"

namespace mpi {


// Exception which is thrown every time a communication fails
struct comm_error : public std::logic_error {
    comm_error(const std::string& msg) :
        std::logic_error(msg) { }
};

class comm;
class endpoint;
class status;
template <class T>
class request;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// comm: is the abstraction of the MPI_Comm class
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class comm{
    MPI_Comm comm_m;
    int comm_size;

    // Check whether MPI_Init has been called
    static void check_init() {
        int flag;
        MPI_Initialized(&flag);
        assert(flag != 0 &&
                "FATAL: MPI environment not initialized (MPI_Init not called)");
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
    comm(MPI_Comm comm): comm_m(comm), comm_size(-1) { }

    // MPI_COMM_WORLD
    static comm world;

    inline size_t rank() {
        check_init();

        int out_rank;
        MPI_Comm_rank(comm_m, &out_rank);
        return out_rank;
    }

    inline size_t size() {
        check_init();

        // Get the size for this communicator
        if(comm_size == -1) {
            MPI_Comm_size(comm_m, &comm_size);
        }
        return comm_size;
    }

    inline bool isMaster() {
        check_init();

        if (rank()==0) { return true; }
        else           { return false; }
    }

    //Converts comm class into an MPI_Comm
    operator MPI_Comm() {return comm_m;}

    inline endpoint operator()( const int& rank_id );

    //Later: Still need to make both send and receive functions
    template <class RawType>
    inline int send (RawType& buffer, int destination, int tag=0);

    template <class RawType>
    inline int recv (RawType& buffer, int source, int tag=0);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Gather functionality for MPP types (types with mpi_type_traits)
//
// RawType1 should be some MPP type T or an MPP container of MPP type T.
// RawType2 should be an MPP container type of the same MPP type T.
// When gathering containers, it is assumed that
//     the receive-buffer-size / comm-size = number of elements to collect
//     from each individual core.
// sendBuffer should be an MPP type that is to be sent.
// recvBuffer should be a container of the same MPP type as what is being sent.
// root is the rank of the core that should receive the message.
// (Containers can be different types, but contained type must be the same.)
// (recvBuffer can be empty on sending-only cores)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template <class RawType1, class RawType2>
    inline int gather (RawType1& sendBuffer, RawType2& recvBuffer, size_t root) {
        msg_impl<RawType1> m(sendBuffer);
        msg_impl<RawType2> r(recvBuffer);
        datatype<RawType1> m_dt(sendBuffer);
        datatype<RawType2> r_dt;
        // gatherv must be used for dynamic types, since each instance is really a different type.
        if (!m.type_static()) {
            std::vector<int> recvCount(size(), recvBuffer.size()/size());
            return gatherv(sendBuffer, recvBuffer, recvCount, root);
        }

        // Check that receive buffer is the correct size.
        if (rank() == root) {
            assert(r.size() == m.size() * size());
            r_dt.reset(recvBuffer);
        }
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
// RawType1 and RawType2 should be MPP container types of some MPP type T.
// sendBuffer should be a container of an MPP type where the entire container
//     is to be sent.
// recvBuffer should be a container of the same MPP type as what is being sent.
// recvCount should be a vector with the Nth element representing how many
//     elements the root should expect to receive from the Nth core.
// root is the rank of the core that should receive the message.
// (Containers can be different types, but contained type must be the same.)
// (recvBuffer and recvCount can be empty on sending-only cores)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template <class RawType1, class RawType2>
    inline int gatherv (RawType1& sendBuffer, RawType2& recvBuffer,
                        const std::vector<int>& recvCount, unsigned int root) {
        msg_impl<RawType1> m(sendBuffer);
        msg_impl<RawType2> r(recvBuffer);
        size_t c_size = this->size();
        int mpi_status;
        std::vector<int>          m_count (c_size);
        std::vector<int>          r_disps (c_size+1);
        m_count[root] = m.size();
        datatype<RawType1> m_type_root(sendBuffer);
        if (m.type_static()) {
            // Check that receive buffer is the correct size.
            if (rank() == root) {
                int totalRecvCount = 0;
                for (size_t i=0; i<c_size; i++) {
                    totalRecvCount += recvCount[i];
                }
                assert(r.size() == totalRecvCount * size());
            }
            datatype<RawType2> r_type_root(recvBuffer);
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
            boost::scoped_array<datatype<RawType2> > r_types(new datatype<RawType2>[c_size]);
            std::vector<int> m_disps (c_size);
            std::vector<int> r_count (c_size);
            int typeDisp=0;
            if (rank()==root) {
                for (size_t i=0; i<c_size; i++) {
                    r_count[i] = recvCount[i] == 0 ? 0 : 1;
                    r_types[i].reset(r.type_range (typeDisp, recvCount[i]+typeDisp));
                    typeDisp  += recvCount[i];
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
// RawType1 and RawType2 should be MPP container types of some MPP type T.
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
    template <class RawType1, class RawType2>
    inline int gatherv (RawType1& sendBuffer, const std::vector<int>& sendIndex,
                        RawType2& recvBuffer, const std::vector<int>& recvCount,
                        const std::vector<int>& recvIndex, unsigned int root) {
        msg_impl<RawType1> m(sendBuffer);
        msg_impl<RawType2> r(recvBuffer);
        size_t c_size = this->size();
        int mpi_status;
        std::vector<int>          m_count (c_size);
        std::vector<int>          r_count (c_size);
        std::vector<int>          m_disps (c_size);
        std::vector<int>          r_disps (c_size);
        boost::scoped_array<datatype<RawType2> > r_types(new datatype<RawType2>[c_size]);
        boost::scoped_array<MPI_Datatype> m_types(new MPI_Datatype[c_size]);
        std::fill(m_types.get(), m_types.get()+c_size, MPI_CHAR);
        m_count[root] = sendIndex.size() == 0 ? 0 : 1;
        assert (check_selection (sendIndex, sendBuffer.size()));
        assert (check_selection (recvIndex, recvBuffer.size()));
        datatype<RawType1> m_type_root(m.type_selection(sendIndex, 0, sendIndex.size()));
        m_types[root] = m_type_root.get();
        if (rank()==root) {
            size_t typeDisp=0;
            for (size_t i=0; i<c_size; i++) {
                r_types[i].reset(r.type_selection(recvIndex, typeDisp, recvCount[i]+typeDisp));
                typeDisp += recvCount[i];
                r_count[i] = recvCount[i] == 0 ? 0 : 1;
            }
        }
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// endpoint: represent the src or dest of an MPI channel. Provides streaming
// operations to send/recv messages (msg<T>) both in a synchronous or asynch
// ronous way.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class endpoint {
    const int           m_rank;  // The rank of this endpoint
    const MPI_Comm&     m_comm;  // The MPI communicator this endpoint
    // belongs to

public:
    endpoint(const int& rank, const MPI_Comm& com):
        m_rank(rank), m_comm(com) { }

    // Send a generic message to this endpoint (synchronously)
    template <class MsgType>
    inline endpoint& operator<<(const msg_impl<MsgType>& m) {
        datatype<MsgType> dt(m.get());
        if (MPI_Send( const_cast<void*>(m.addr()), m.size(), dt.get(),
                        m_rank, m.tag(), m_comm
                    ) == MPI_SUCCESS ) {
            return *this;
        }
        std::ostringstream ss;
        ss << "ERROR in MPI rank '" << comm::world.rank()
            << "': Failed to send message to destination rank '"
            << m_rank << "'";
        throw comm_error( ss.str() );
    }

    // Send a generic raw type to this endpoint (synchronously)
    template <class RawType>
    inline endpoint& operator<<(RawType& m) {
        return operator<<(msg_impl<RawType>(m));
    }

    // Receive from this endpoint (synchronously)
    template <class RawType>
    inline status operator>>(RawType& m);

    template <class MsgType>
    inline status operator>>(const msg_impl<MsgType>& m);

    // Receive from this endpoint (asynchronously)
    template <class MsgType>
    inline request<MsgType> operator>(const msg_impl<MsgType>& m) {
        MPI_Request req;
        datatype<MsgType> dt(m.get());
        if(MPI_Irecv( const_cast<void*>(m.addr()), m.size(), dt,
                        m_rank, m.tag(), m_comm, &req
                    ) != MPI_SUCCESS ) {
            std::ostringstream ss;
            ss << "ERROR in MPI rank '" << comm::world.rank()
                << "': Failed to receive message from destination rank '"
                << m_rank << "'";
            throw comm_error( ss.str() );
        }
        return request<MsgType>(m_comm, req, m);
    }

    // Receive from this endpoint (asynchronously
    template <class RawType>
    inline request<RawType> operator>(RawType& m) {
        return operator>( msg_impl<RawType>( m ) );
    }

    // Returns the rank of this endpoint
    inline const int& rank() const { return m_rank; }
};

endpoint comm::operator()(const int& rank_id) {
    return endpoint(rank_id, comm_m);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// status: keeps the status info for received messages.
// containes the sender of the message, tag and size of
// the message
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class status{
    const MPI_Comm&      m_comm;
    const MPI_Status     m_status;
public:
    status(const MPI_Comm& com, const MPI_Status& s):
            m_comm(com), m_status(s)
         { }

    inline endpoint source(){
        return endpoint(m_status.MPI_SOURCE, m_comm);
    }

    inline int tag(){
        return m_status.MPI_TAG;
    }

    inline int error(){
        return m_status.MPI_ERROR;
    }
};

template <class RawType>
inline status endpoint::operator>>(RawType& m) {
    return operator>>( msg_impl<RawType>( m ) );
}

template <class MsgType>
inline status endpoint::operator>>(const msg_impl<MsgType>& m) {
    MPI_Status s;
    datatype<MsgType> dt(m.get());
    if (MPI_Recv( const_cast<void*>(m.addr()), m.size(), dt.get(),
                m_rank, m.tag(), m_comm, &s
                ) == MPI_SUCCESS ) {
        return status(m_comm, s);
    }
    std::ostringstream ss;
    ss << "ERROR in MPI rank '" << comm::world.rank()
       << "': Failed to receive message from destination rank '"
       << m_rank << "'";
    throw comm_error( ss.str() );
}

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
