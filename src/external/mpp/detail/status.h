
#pragma once

#include "decls.h"

namespace mpi
{

class status{
    const MPI_Comm&      m_comm;
    const MPI_Status     m_status;
public:
    status(const MPI_Comm& com, const MPI_Status& s):
            m_comm(com), m_status(s)
         { }

//    inline endpoint source(){  //TODO: add back in (removed because of cyclic header dependencies)
//        return endpoint(m_status.MPI_SOURCE, m_comm);
//    }

    inline int tag(){
        return m_status.MPI_TAG;
    }

    inline int error(){
        return m_status.MPI_ERROR;
    }
};

} // end mpi namespace
