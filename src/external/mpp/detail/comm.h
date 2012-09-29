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
#include <numeric>
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
                          << "' but tried to access index '" << s << "'\n";
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

    comm(const comm &other):
        m_initialized(other.m_initialized),
        m_comm_size(other.m_comm_size),
        m_rank(other.m_rank)
    { MPI_Comm_dup(other.m_comm, &m_comm); }

    ~comm() {
        if(m_comm!=world.m_comm)
            MPI_Comm_free(&m_comm);
    }

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
    inline void send (RawType& buffer, int destination, int tag=0);

    template <class RawType>
    inline void recv (RawType& buffer, int source, int tag=0);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Gather functionality for MPP types (types with mpi_type_traits)
//
// When gathering containers, it is assumed that
//     the receive-buffer-size / comm-size = number of elements to collect
//     from each individual core.
// root is the rank of the core that should receive the message.
// Iterator pairs describing a range need to be of the same type.
// Iterators of different pairs can be of different types and also point
// to different types, as long as the underlying primitive is the same.
// E.g. it is OK to send a std::list<struct{int}>::iterator and receive
// into a int pointer.
// recvBuffer can be empty on sending-only cores
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template <class SendIterator, class RecvIterator>
    inline void gather (SendIterator sendBegin, SendIterator sendEnd,
            RecvIterator recvBegin, RecvIterator recvEnd, int root) {
#ifndef NDEBUG
        // Check that receive buffer is the correct size.
        if (rank() == root) {
        //not ok for e.g. VectorObjToVector. Can we check whether size matches?
//            assert(std::distance(recvBegin,recvEnd) == nElem * this->size());
        }
#endif
        if (sendBegin==sendEnd) {
            return; //nothing to do
        }
        // gatherv must be used for non-contiguous recv containers, since each instance is really a different type.
        if (!is_contiguous<RecvIterator>::value) {
            std::vector<int> recvCount;
            if (rank() == root) {
                recvCount.resize(size());
                int nElem = std::distance(recvBegin,recvEnd);
                assert(nElem%this->size()==0);
                std::fill(recvCount.begin(), recvCount.end(), nElem/this->size());
            }
            return gatherv(sendBegin, sendEnd, recvBegin, recvEnd,
                    recvCount.begin(), recvCount.end(), root);
        }
        data_layout recv_dt;
        if (rank() == root) {
            recv_dt = get_layout(recvBegin, recvEnd);
        }
        data_layout send_dt = get_layout(sendBegin, sendEnd);
        MPI_Gather (
            send_dt.get_ptr(), send_dt.get_size(), send_dt.get(),
            recv_dt.get_ptr(), recv_dt.get_size()/this->size(), recv_dt.get(), root, *this);
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Gatherv functionality for MPP types (types with mpi_type_traits)
//
// [sendBegin, sendEnd) range is sent.
// [recvBegin, recvEnd) range is received into.
// [recvCountBegin, recvCountEnd) range gives the number of elemets received
//       from each rank. The length of the range has to be comm.size()
// The underlying types send/recv ranges need to be compatible.
// root is the rank of the core that should receive the message.
// recvBuffer and recvCount can be empty on sending-only cores
// CountIterator has to be contiguous and its value_type has to be int
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template <class SendIterator, class RecvIterator, class CountIterator>
    inline void gatherv (SendIterator sendBegin, SendIterator sendEnd,
            RecvIterator recvBegin, RecvIterator recvEnd,
            CountIterator recvCountBegin, CountIterator recvCountEnd,
            int root) {
        int comm_size = this->size();
        std::vector<int>          recv_disps (comm_size+1);
        //datatype for whole send range (for send no is_contiguous check is done)
        if (is_contiguous<RecvIterator>::value) {
            //datatype for single element in recv range (OK because is_contiguous)
            data_layout recv_type;
            std::vector<int> real_recv_count(comm_size); //corrected for recv_type.get_size()
            if (rank() == root) {
                // Check that receive buffer is the correct size.
                assert(std::distance(recvBegin,recvEnd) ==
                        std::accumulate(recvCountBegin, recvCountEnd, 0));
//                assert(std::distance(sendBegin,sendEnd) == *recvCountBegin);
                assert(recvCountEnd-recvCountBegin == comm_size);
                assert(is_contiguous<CountIterator>::value);

                if (recvBegin!=recvEnd) {
                    recv_type = get_layout(*recvBegin);
                    for(int i=0; i<comm_size; i++) {
                        real_recv_count[i] = recvCountBegin[i]*recv_type.get_size();
                        recv_disps[i+1] = recv_disps[i] + real_recv_count[i];
                    }
                }
            }
            data_layout send_type = get_layout(sendBegin, sendEnd);
            MPI_Gatherv (
                send_type.get_ptr(), send_type.get_size(), send_type.get(),
                recv_type.get_ptr(), &real_recv_count.front(),
                &recv_disps.front(), recv_type.get(), root, *this);
        } else {
            std::vector<MPI_Datatype> send_types(comm_size, MPI_CHAR);
            std::vector<int> send_count (comm_size);
            std::vector<int> send_disps (comm_size);
            data_layout send_type_root = get_layout(sendBegin, sendEnd);
            send_types[root] = send_type_root.get();
            send_count[root] = send_type_root.get_size();

            std::vector<MPI_Datatype> recv_types(comm_size);
            boost::scoped_array<data_layout> recv_ranges(new data_layout[comm_size]);
            std::vector<int> recv_count (comm_size);

            if (rank()==root) {
                RecvIterator typeDisp = recvBegin;
                for (int i=0; i<comm_size; i++) {
                    RecvIterator typeEnd = typeDisp;
                    std::advance(typeEnd, recvCountBegin[i]);
                    //needs to use dt with absolute address
                    recv_ranges[i] = get_layout(typeDisp, typeEnd, true);
                    recv_types[i] = recv_ranges[i].get();
                    recv_count[i] = recv_ranges[i].get_size();
                    typeDisp = typeEnd;
                }
            }
            //TODO: replace with own GatherW. OpenMPI(1.6) & MPICH2(1.5) anyhow only do send/irecv even for gatherv
            //(maybe fallback to alltoallw for intercomm)

            // m_disps and r_disps have all zero values. r_count has all zero values for non-root and for root has all values
            // one or zero (if receive count is zero). m_count is zero for non-root and is the size of the send buffer for root.
            MPI_Alltoallw(
                send_type_root.get_ptr(), &send_count.front(), &send_disps.front(), &send_types.front(),
                MPI_BOTTOM, &recv_count.front(), &recv_disps.front(), &recv_types.front(), *this);
        }
    }

    template <class SendIterator, class RecvIterator>
    inline void reduce (SendIterator sendBegin, SendIterator sendEnd,
            RecvIterator recvBegin, RecvIterator recvEnd, MPI_Op op, int root) {
        int nElem = std::distance(sendBegin,sendEnd);
        if (nElem==0) {
            return; //nothing to do
        }
        //currently only supports contiguous or size of 1 because MPI doesn't has the standard operations
        //defined for derived datatypes (not even MPI_Type_contigouus)
        assert(is_contiguous<SendIterator>::value || nElem==1);
        bool inplace;
        // Check that receive buffer is the correct size.
        if (rank() == root) {
            int nElemRecv = std::distance(recvBegin,recvEnd);
            inplace = (nElemRecv==0);
            assert(is_contiguous<RecvIterator>::value || nElemRecv==1);
            assert( nElemRecv == nElem || inplace);
        }
        data_layout send_dt = get_layout(sendBegin,sendEnd);
#ifndef NDEBUG
        if (rank() == root && !inplace) {
            data_layout recv_dt = get_layout(*recvBegin);
            //because above is tested that type is static the datatype it is correct to test
            //for equal (and not necessary to test for identical)
            assert(recv_dt.get()==send_dt.get());
        }
#endif
        void *sendaddr, *recvaddr = MPI_BOTTOM;
        sendaddr = send_dt.get_ptr();
        if (rank() == root) {
            if (inplace) {
                recvaddr = sendaddr;
                sendaddr = MPI_IN_PLACE;
            } else {
                recvaddr = &*recvBegin; //for static safe
            }
        }
        MPI_Reduce (sendaddr, recvaddr, send_dt.get_size(), send_dt.get(), op, root, *this);
    }

    //inplace version
    template <class Iterator>
        inline void reduce (Iterator begin, Iterator end, MPI_Op op, int root) {
            reduce(begin,end,(char*)NULL,(char*)NULL,op,root); //type of NULL has to be some pointer
    }
};

inline void init(int* argc = 0, char** argv[] = NULL){
    MPI_Init(argc, argv);
    comm::world = comm(MPI_COMM_WORLD);
}
inline void finalize() {
    cached_types.clear();
    MPI_Finalize();
}
const int any = MPI_ANY_SOURCE;

} // end mpi namespace
