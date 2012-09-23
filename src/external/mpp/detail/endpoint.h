
#pragma once

#include "decls.h"

namespace mpi
{
//TODO: Add (a non operator version of) Send/Recv accepting iterators.
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
        if (MPI_Send( const_cast<void*>(m.addr()), 1, dt.get(),
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
    inline request<MsgType> operator>(const msg_impl<MsgType>& m) {
        MPI_Request req;
        datatype<MsgType> dt(m.get());
        if(MPI_Irecv( m.addr(), 1, dt.get(),
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

} // end mpi namespace

#include "status.h"

namespace mpi
{

endpoint comm::operator()(const int& rank_id) {
    return endpoint(rank_id, m_comm);
}


template <class RawType>
inline status endpoint::operator>>(RawType& m) {
    return operator>>( msg_impl<RawType>( m ) );
}

template <class MsgType>
inline status endpoint::operator>>(const msg_impl<MsgType>& m) {
    MPI_Status s;
    datatype<MsgType> dt(m.get());
    if (MPI_Recv( const_cast<void*>(m.addr()), 1, dt.get(),
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


} // end mpi namespace
