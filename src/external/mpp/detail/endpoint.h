
#pragma once

#include "decls.h"

namespace mpi
{
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
        MPI_Send( const_cast<void*>(m.addr()), m.size(), m.type(),
                        m_rank, m.tag(), m_comm);
        return *this;
    }

    // Send a generic raw type to this endpoint (synchronously)
    template <class RawType>
    inline endpoint& operator<<(const RawType& m) {
        return operator<<(msg_impl<const RawType>(m));
    }

    // Receive from this endpoint (synchronously)
    template <class RawType>
    inline status operator>>(RawType& m);

    template <class MsgType>
    inline status operator>>(const msg_impl<MsgType>& m);

    // Receive from this endpoint (asynchronously)
    template <class MsgType>
    inline request operator>(const msg_impl<MsgType>& m);

    // Receive from this endpoint (asynchronously
    template <class RawType>
    inline request operator>(RawType& m);

    // Returns the rank of this endpoint
    inline const int& rank() const { return m_rank; }
};

} // end mpi namespace

#include "status.h"
#include "request.h"
#include "comm.h"

namespace mpi
{

template <class MsgType>
inline request endpoint::operator>(const msg_impl<MsgType>& m) {
    MPI_Request req;
    MPI_Irecv( m.addr(), m.size(), m.type(),
               m_rank, m.tag(), m_comm, &req);
    return request(m_comm, req);
}

template <class RawType>
inline request endpoint::operator>(RawType& m) {
    return operator>( msg_impl<RawType>( m ) );
}

endpoint comm::operator()(const int& rank_id) {
    return endpoint(rank_id, m_comm);
}


template <class RawType>
inline status endpoint::operator>>(RawType& m) {
    return operator>>( msg_impl<RawType>( m ) );
}

//TODO: should call comm::recv
template <class MsgType>
inline status endpoint::operator>>(const msg_impl<MsgType>& m) {
    MPI_Status s;
    MPI_Recv( const_cast<void*>(m.addr()), m.size(), m.type(),
                m_rank, m.tag(), m_comm, &s);
    return status(m_comm, s);
}


} // end mpi namespace
