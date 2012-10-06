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

#include <limits.h>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "message.h"
#include "decls.h"
const int COLLECTIVE_TAG = INT_MAX-1337;

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

    /*return remote size if inter-communicator otherwise returns size */
    inline int remote_size() const {
        if (is_inter_comm()) {
            int size;
            MPI_Comm_remote_size(m_comm,&size);
            return size;
        } else {
            return size();
        }
    }

    inline bool is_inter_comm() const {
        int flag;
        MPI_Comm_test_inter(m_comm, &flag);
        return flag;
    }

    inline bool isMaster() {
        check_init();

        if (rank()==0) { return true; }
        else           { return false; }
    }

    //Converts comm class into an MPI_Comm
    operator MPI_Comm() {return m_comm;}

    inline endpoint operator()( const int& rank_id );

    //Wrapper around send (itr, itr,....)
    template <class RawType>
    inline void send (RawType& buffer, int destination, int tag=0) {
        send (&buffer, (&buffer)+1, destination, tag);
    }

    template <class startIterator, class endIterator>
    inline void send (startIterator& sendBegin, endIterator& sendEnd, int destination, int tag=0) {
        data_layout send_type = get_range_layout(sendBegin, sendEnd);
        MPI_Send (send_type.get_ptr(), send_type.get_size(), send_type.get(), destination, tag, *this);
    }

    //Wrapper around recv (itr, itr,....)
    template <class RawType>
    inline status recv (RawType& buffer, int source, int tag=0) {
        return recv (&buffer, (&buffer)+1, source, tag);
    }

    template <class startIterator, class endIterator>
    inline status recv (startIterator& recvBegin, endIterator& recvEnd, int source, int tag=0) {
        MPI_Status s;
        data_layout recv_type = get_range_layout(recvBegin, recvEnd);
        MPI_Recv (recv_type.get_ptr(), recv_type.get_size(), recv_type.get(), source, tag, *this, &s);
        return status(*this, s);
    }

    //Wrapper around irecv (itr, itr,....)
    template <class RawType>
    inline request irecv (RawType& buffer, int source, int tag=0) {
        return irecv (&buffer, (&buffer)+1, source, tag);
    }

    template <class startIterator, class endIterator>
    inline request irecv (startIterator& recvBegin, endIterator& recvEnd, int source, int tag=0) {
        MPI_Request r;
        data_layout recv_type = get_range_layout(recvBegin, recvEnd);
        MPI_Irecv (recv_type.get_ptr(), recv_type.get_size(), recv_type.get(), source, tag, *this, &r);
        return request(*this, r);
    }
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
                recvCount.resize(this->size());
                int nElem = std::distance(recvBegin,recvEnd);
                assert(nElem%this->size()==0);
                std::fill(recvCount.begin(), recvCount.end(), nElem/this->size());
            }
            return gatherv(sendBegin, sendEnd, recvBegin, recvEnd,
                    recvCount.begin(), recvCount.end(), root);
        }
        data_layout recv_dl;
        if (rank() == root) {
            recv_dl = get_range_layout(recvBegin, recvEnd);
        }
        data_layout send_dl = get_range_layout(sendBegin, sendEnd);
        MPI_Gather (
            send_dl.get_ptr(), send_dl.get_size(), send_dl.get(),
            recv_dl.get_ptr(), recv_dl.get_size()/this->size(), recv_dl.get(), root, *this);
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
        if (rank() == is_inter_comm()?MPI_ROOT:root) {
            // Check that receive buffer is the correct size.
            assert(std::distance(recvBegin,recvEnd) ==
                    std::accumulate(recvCountBegin, recvCountEnd, 0));
    //                assert(std::distance(sendBegin,sendEnd) == *recvCountBegin);
            assert(recvCountEnd-recvCountBegin == remote_size());

            std::vector<request> requests;
            RecvIterator sectionBegin = recvBegin, sectionEnd = recvBegin;
            for (int i=0; i<remote_size(); i++) {
                std::advance(sectionEnd, recvCountBegin[i]);
                request r = irecv (sectionBegin, sectionEnd, i, COLLECTIVE_TAG);
                requests.push_back (move(r));
                sectionBegin = sectionEnd;
            }
            if (!is_inter_comm()) {
                send (sendBegin, sendEnd, root, COLLECTIVE_TAG);
            }
            wait (requests.begin(), requests.end());
        } else if (root != MPI_PROC_NULL) { /* non-root nodes, and in the intercomm. case, non-root nodes on remote side */
            send (sendBegin, sendEnd, root, COLLECTIVE_TAG);
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
        data_layout send_dl = get_range_layout(sendBegin,sendEnd);
#ifndef NDEBUG
        if (rank() == root && !inplace) {
            data_layout recv_dl = get_type_layout(*recvBegin);
            //because above is tested that type is static the datatype it is correct to test
            //for equal (and not necessary to test for identical)
            assert(recv_dl.get()==send_dl.get());
        }
#endif
        void *sendaddr, *recvaddr = MPI_BOTTOM;
        sendaddr = send_dl.get_ptr();
        if (rank() == root) {
            if (inplace) {
                recvaddr = sendaddr;
                sendaddr = MPI_IN_PLACE;
            } else {
                recvaddr = &*recvBegin; //for static safe
            }
        }
        MPI_Reduce (sendaddr, recvaddr, send_dl.get_size(), send_dl.get(), op, root, *this);
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
    cached_layouts.clear();
    MPI_Finalize();
}
const int any = MPI_ANY_SOURCE;

} // end mpi namespace
