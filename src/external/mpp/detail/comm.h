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
    inline int send (RawType& buffer, int destination, int tag=0);

    template <class RawType>
    inline int recv (RawType& buffer, int source, int tag=0);

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
    inline int gather (SendIterator sendBegin, SendIterator sendEnd,
            RecvIterator recvBegin, RecvIterator recvEnd, int root) {
        // gatherv must be used for non-contiguous recv containers, since each instance is really a different type.
        int nElem = std::distance(sendBegin,sendEnd);
#ifndef NDEBUG
        // Check that receive buffer is the correct size.
        if (rank() == root) {
            assert(std::distance(recvBegin,recvEnd) == nElem * this->size());
        }
#endif
        if (!is_contiguous<RecvIterator>::value) {
            std::vector<int> recvCount;
            if (rank() == root) {
                recvCount.resize(size());
                int nElem = std::distance(recvBegin,recvEnd)/size();
                std::fill(recvCount.begin(), recvCount.end(), nElem);
            }
            return gatherv(sendBegin, sendEnd, recvBegin, recvEnd,
                    recvCount.begin(), recvCount.end(), root);
        }

        if (sendBegin==sendEnd) {
            return MPI_SUCCESS; //nothing to do
        }
        typedef typename std::iterator_traits<SendIterator>::value_type SendType;
        typedef typename std::iterator_traits<RecvIterator>::value_type RecvType;
        datatype<RecvType> recv_dt;
        void* recv_addr = MPI_BOTTOM;
        if (rank() == root) {
            recv_dt.reset(*recvBegin);
            recv_addr = MPP_Get_ptr(*recvBegin);
        }

        int ret;
        if (nElem == 1 || is_contiguous<SendIterator>::value ) {
            ret = MPI_Gather (MPP_Get_ptr(*sendBegin), nElem,
                datatype<SendType>(*sendBegin).get(),
                recv_addr, nElem, recv_dt.get(), root, *this);
        } else {
            ret = MPI_Gather (MPI_BOTTOM, nElem>0?1:0,
                datatype_range(sendBegin, sendEnd).get(),
                recv_addr, nElem, recv_dt.get(), root, *this);
        }
        if (ret != MPI_SUCCESS) {
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
    inline int gatherv (SendIterator sendBegin, SendIterator sendEnd,
            RecvIterator recvBegin, RecvIterator recvEnd,
            CountIterator recvCountBegin, CountIterator recvCountEnd,
            int root) {
        typedef typename std::iterator_traits<SendIterator>::value_type SendType;
        typedef typename std::iterator_traits<RecvIterator>::value_type RecvType;
        int comm_size = this->size();
        int mpi_status;
        std::vector<int>          resv_disps (comm_size);
        //datatype for whole send range (for send no is_contiguous check is done)
        if (is_contiguous<RecvIterator>::value) {
            //datatype for single element in recv range (OK because is_contiguous)
            datatype<RecvType> recv_type_root;
            void* recv_addr = MPI_BOTTOM;
            if (rank() == root) {
                // Check that receive buffer is the correct size.
                assert(std::distance(recvBegin,recvEnd) ==
                        std::accumulate(recvCountBegin, recvCountEnd, 0));
                assert(std::distance(sendBegin,sendEnd) == *recvCountBegin);
                assert(recvCountEnd-recvCountBegin == comm_size);
                assert(is_contiguous<CountIterator>::value);

                if (recvBegin!=recvEnd) {
                    recv_type_root.reset(*recvBegin);
                    recv_addr = MPP_Get_ptr(*recvBegin);
                } //no else required - using initialization value
                //set resv_disps[1,end) to partial sum of [recvCountBegin,recvCountEnd-1)
                std::partial_sum(recvCountBegin, --recvCountEnd, ++(resv_disps.begin()));
            }
            //TODO: It would be OK to always use datatype_range but it creates a hindex dt
            //which is not necessary if range is contiguous or we only send 1 element.
            //This optimization should be done for every send operation which benefits.
            //Thus it should be abstracted out so that this if/else code currently here
            //and in Gather isn't required for each operation. And then it should be added
            //to the Alltoallw branch here in Gatherv.
            SendIterator sendBeginPlus1 = sendBegin;
            ++sendBeginPlus1;
            if (sendBeginPlus1==sendEnd || is_contiguous<SendIterator>::value ) {
                mpi_status = MPI_Gatherv (MPP_Get_ptr(*sendBegin),
                    std::distance(sendBegin,sendEnd),
                    datatype<SendType>(*sendBegin).get(),
                    recv_addr, const_cast<int*>(&*recvCountBegin),
                    &resv_disps.front(), recv_type_root.get(), root, *this);
            } else {
                mpi_status = MPI_Gatherv (MPI_BOTTOM, (sendBegin==sendEnd)?0:1,
                    datatype_range(sendBegin, sendEnd).get(),
                    recv_addr, const_cast<int*>(&*recvCountBegin),
                    &resv_disps.front(), recv_type_root.get(), root, *this);
            }
        } else {
            boost::scoped_array<MPI_Datatype> send_types(new MPI_Datatype[comm_size]);
            std::fill(send_types.get(), send_types.get()+comm_size, MPI_CHAR); //Using 0 chars for empty messages
            datatype_range send_type_root(sendBegin, sendEnd);
            send_types[root] = send_type_root.get();
            boost::scoped_array<datatype_range> recv_types(new datatype_range[comm_size]);
            std::vector<int> send_disps (comm_size);
            std::vector<int> recv_count (comm_size);
            std::vector<int> send_count (comm_size);
            send_count[root] = (sendBegin==sendEnd)?0:1;
            if (rank()==root) {
                RecvIterator typeDisp = recvBegin;
                for (int i=0; i<comm_size; i++) {
                    recv_count[i] = recvCountBegin[i] == 0 ? 0 : 1;
                    RecvIterator typeEnd = typeDisp;
                    std::advance(typeEnd, recvCountBegin[i]);
                    recv_types[i].reset(typeDisp, typeEnd);
                    typeDisp = typeEnd;
                }
            }
            // m_disps and r_disps have all zero values. r_count has all zero values for non-root and for root has all values
            // one or zero (if receive count is zero). m_count is zero for non-root and is the size of the send buffer for root.
            mpi_status = MPI_Alltoallw(
                    MPI_BOTTOM, &send_count.front(), &send_disps.front(), send_types.get(),
                    MPI_BOTTOM, &recv_count.front(), &resv_disps.front(),
                    const_cast<MPI_Datatype*>(&recv_types.get()->get()), *this);
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

    template <class SendIterator, class RecvIterator>
    inline int reduce (SendIterator sendBegin, SendIterator sendEnd,
            RecvIterator recvBegin, RecvIterator recvEnd, MPI_Op op, int root) {
        //currently only supports contiguous because MPI doesn't has the standard operations
        //defined for derived datatypes (not even MPI_Type_contigouus)
        assert(is_contiguous<RecvIterator>::value && is_contiguous<SendIterator>::value);
        int nElem = std::distance(sendBegin,sendEnd);
        bool inplace;
        // Check that receive buffer is the correct size.
        if (rank() == root) {
            int nElemRecv = std::distance(recvBegin,recvEnd);
            inplace = (nElemRecv==0);
            assert( nElemRecv == nElem || inplace);
        }
        if (nElem==0) {
            return MPI_SUCCESS; //nothing to do
        }
        typedef typename std::iterator_traits<SendIterator>::value_type SendType;
        typedef typename std::iterator_traits<RecvIterator>::value_type RecvType;
        datatype<SendType> send_dt(*sendBegin);
#ifndef NDEBUG
        if (rank() == root && !inplace) {
            datatype<RecvType> recv_dt(*recvBegin);
            //because above is tested that type is static the datatype it is correct to test
            //for equal (and not necessary to test for identical)
            assert(recv_dt.get()==send_dt.get());
        }
#endif
        void *sendaddr, *recvaddr = MPI_BOTTOM;
        sendaddr = (void*)MPP_Get_ptr(*sendBegin);
        if (rank() == root) {
            if (inplace) {
                recvaddr = sendaddr;
                sendaddr = MPI_IN_PLACE;
            } else {
                recvaddr = MPP_Get_ptr(*recvBegin);
            }
        }
        if (MPI_Reduce (sendaddr, recvaddr, nElem, send_dt.get(), op, root, *this)
                        != MPI_SUCCESS) {
            std::ostringstream ss;
            ss << "ERROR in MPI rank '" << comm::world.rank()
               << "': Failed to gather message to destination rank '"
               << root << "'";
            throw comm_error( ss.str() );
        }
        return MPI_SUCCESS;
    }

    //inplace version
    template <class Iterator>
        inline int reduce (Iterator begin, Iterator end, MPI_Op op, int root) {
            return reduce(begin,end,(char*)NULL,(char*)NULL,op,root); //type of NULL has to be some pointer
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
