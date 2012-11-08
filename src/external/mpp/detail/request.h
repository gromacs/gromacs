#pragma once

#include "decls.h"

namespace mpi {
class request {
    public:
        request(MPI_Comm comm, MPI_Request req) :
            storage_(new request_storage(comm, req)) {}
        status wait();
    private:
        struct request_storage {
            request_storage(MPI_Comm comm, MPI_Request req) :
                req_(req), comm_(comm) {};
            ~request_storage() {
                if(req_!=MPI_REQUEST_NULL) {
                    MPI_Request_free(&req_);
                }
            }
            MPI_Request req_;
            MPI_Comm comm_;
        };
        mpp_unique_ptr<request_storage>::type storage_;
};

#ifndef MPP_CXX11_RVALREF
static inline request &move(request &r) { return r; } //No-Op
#endif

//wait for all requests in range [begin,end)
template <class Iter>
std::vector<status> wait(Iter begin, Iter end);

} // end mpi namespace

#include "status.h"

namespace mpi
{

inline  status request::wait() {
    MPI_Status s;
    MPI_Wait(&(storage_->req_), &s);
    return status(storage_->comm_, s);
}

template <class Iter>
std::vector<status> wait(Iter begin, Iter end) {
    std::vector<status> statuses;
    statuses.reserve(std::distance(begin, end));
    for (Iter i=begin; i!=end; i++) {
        statuses.push_back(i->wait());
    }
    return statuses;
}
} //end namespace mpi
