#ifndef TASK_2D_H
#define TASK_2D_H

// TODO #include <typedefs.h>
#include "progress.h"
#include "utils.hpp"

namespace FCA {

class Task2D {
    static const int __TAG       = 42;
    static const int REQUEST_TAG = 44;

    const utils::mpi& com;
    double* mat;
    int N, M;
    int count_row;
    int count_col;
    int actual_row;
    int actual_col;
    int bSym;
    double result;
    const std::pair<int,int>* pairs; /* if set, only compute these pairs */
    int pair_count;
    real* to_compute; /* compute value only if entry in to_compute larger than "compute_threshold" */
    real compute_threshold;
    Progress prog;

#ifdef MPICC
    static void _receive(void* data, int nr, MPI_Datatype type, int address, int tag, MPI_Comm comm, MPI_Status* stat) {
        AssertMpi(MPI_Recv(data, nr, type, address, tag, comm, stat));
    }

    static void _send(void* dat, int nr, MPI_Datatype type, int address, int tag, MPI_Comm comm) {
        /*  SLOG(send_log() << "send_to_" << address << "_from_" << my_rank << std::endl)*/
        AssertMpi(MPI_Ssend(dat, nr, type, address, tag, comm));
        /* SLOG(send_log() << "...done" << "send_to_" << address << "_from_" << my_rank << std::endl) */
    }

    static void mpi_send_it(const utils::mpi& /*mpiInfo*/, int address, int dat, int tag) {
        _send(&dat, 1, MPI_INT, address, tag, MPI_COMM_WORLD);
    }

    static void mpi_send_i(const utils::mpi& mpiInfo, int address, int dat) {
        mpi_send_it(mpiInfo, address, dat, __TAG);
    }

    static void mpi_send_d(const utils::mpi& /*mpiInfo*/, int address, double dat) {
        _send(&dat, 1, MPI_DOUBLE, address, __TAG, MPI_COMM_WORLD);
    }

    static void mpi_receive_d(const utils::mpi& /*mpiInfo*/, int address, double* data) {
        MPI_Status stat;
        _receive(data, 1, MPI_DOUBLE, address, __TAG, MPI_COMM_WORLD, &stat);
    }

    static void mpi_receive_i(const utils::mpi& /*mpiInfo*/, int address, int* data) {
        MPI_Status stat;
        _receive(data, 1, MPI_INT, address, __TAG, MPI_COMM_WORLD, &stat);
    }

    static int mpi_receive_request(const utils::mpi& /*mpiInfo*/) {
        int data;
        MPI_Status stat;
        _receive(&data, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &stat);
        return data;
    }

    static void mpi_send_request(const utils::mpi& mpiInfo) {
        mpi_send_it(mpiInfo, 0, mpiInfo.getMyRank(), REQUEST_TAG);
    }
#endif

public:

#ifdef MPICC
    int receive_result() {
        int slave;
        double data;
        int i, j;
        slave = mpi_receive_request(this->com);
        mpi_receive_d(this->com, slave, &data);
        mpi_receive_i(this->com, slave, &i);
        mpi_receive_i(this->com, slave, &j);
        /* fprintf(stderr,"receive from node %d for cell (%d,%d) %g\n",slave,i,j,data); */
        /* receives one bogus-result for cell (-1,-1) in beginning from every node, do not store... */
        if(i >= 0) {
            this->mat[i * this->M + j] = data;
            if(this->bSym) {
                this->mat[j * this->N + i] = data;
            }
        }
        return slave;
    }

    void send_result() {
        mpi_send_request(this->com);
        mpi_send_d(this->com, 0, this->result);
        mpi_send_i(this->com, 0, this->actual_row);
        mpi_send_i(this->com, 0, this->actual_col);
    }
#endif

    Task2D(const utils::mpi& inCom, double* mat, const int M, const int N, const bool bSym)
    : com(inCom) {
        this->count_col  = 0;
        this->count_row  = 0;
        this->actual_row = -1;
        this->actual_col = -1;
        this->M          = M;
        this->N          = N;
        this->bSym       = bSym;
        this->mat        = mat;
        this->to_compute = nullptr;
        if(bSym && M != N) {
            fprintf(stderr, "init requesting a symmetric matrix with different dimensions %d x %d\n", M, N);
            exit(1);
            /* replaced this due to compatibility issues gmx3.3 vs gmx3.2
               gmx_fatal(FARGS,"init requesting a symmetric matrix with different dimensions %d x %d\n",M,N); */
        }
        if(com.getMyRank() == MASTER_RANK) {
            this->prog = Progress("", (bSym ? N * (N + 1) / 2 : N * this->N),
                                  100,                   /* parallel tasks are very slow ... get more ticks to view */
                                  -com.getNumNodes() + 1 /* in the beginning there is a start tick for all nodes this is no real progress. */
            );
        }
        /* for sym: N*(N+1)/2 since diagonal elements are ticked, too */
        this->pairs      = nullptr;
        this->pair_count = -1;
    }

    Task2D(const utils::mpi& com, double* mat, const int M, const int N, const bool bSym, const std::pair<int,int> pairs[], const int npairs)
            : Task2D(com, mat, M, N, bSym) {
        this->pairs      = pairs;
        this->pair_count = npairs;
    }

    virtual ~Task2D() {
    }

    Task2D(const Task2D&) = delete;
    Task2D& operator=(const Task2D&) = delete;

    void set_to_compute(real* to_compute, real threshold) {
        this->to_compute        = to_compute;
        this->compute_threshold = threshold;
    }

    int finished() {
        if(this->pairs) {
            return this->count_col >= this->pair_count;
        } else {
            return this->count_row >= this->M;
        }
    }

    int newcol() {
        if(this->pairs) {
            // fprintf(stderr,"NEWCOL   count_col %d (%d)   count_row %d \n",this->count_col,this->pair_count,this->count_row);
            if(!finished()) {
                while(this->pairs[this->count_col].first <= 0) {
                    this->count_col++;
                    if(this->count_col >= this->pair_count) {
                        return -1;
                    }
                }
                this->count_row = this->pairs[this->count_col].first;
                // fprintf(stderr,"do move %d (%d %d)\n",
                //        count_col,
                //        this->pairs[this->count_col].first,this->pairs[this->count_col].second);
                return this->pairs[(this->count_col)++].second;
            } else {
                return -1;
            }
        } else {
            if(((this->count_col > this->count_row) && this->bSym) || (this->count_col >= this->N)) {
                this->count_col = 0;
                ++(this->count_row);
            }
            if(!finished()) {
                if(this->to_compute) {
                    if(this->to_compute[this->count_row + this->count_col * this->M] > this->compute_threshold) {
                        this->to_compute[this->count_row + this->count_col * this->M] = 0;
                        return (this->count_col)++;
                    } else {
                        ++this->count_col;
                        return newcol();
                    }
                } else { //that is !to_compute
                    return (this->count_col)++;
                }
            } else {
                return -1;
            }
        }
        return EXIT_FAILURE;
    }


    void push_result(double res) {
        this->result = res;
        if(this->com.getNumNodes() == 1) {
            /* no MPICC or only single node */
            this->mat[this->actual_row * this->M + this->actual_col] = res;
            if(this->bSym) {
                this->mat[this->actual_col * this->N + this->actual_row] = res;
            }
        }
    }

    int next_index(int* i, int* j) {
        if(this->com.getNumNodes() > 1) {
#ifdef MPICC
            if(this->com.getMyRank() == MASTER_RANK) {
                /*one job is for the master himself */
                int col = newcol();
                while(col >= 0) {
                    int slave = receive_result();
                    //       	progress_tick(&this->prog); /* tick before the job is done */
                    mpi_send_i(this->com, slave, col);
                    mpi_send_i(this->com, slave, this->count_row);
                    col = newcol();
                }
                for(int slave = 1; slave < this->com.getNumNodes(); slave++) {
                    int nex_slave = receive_result();
                    //	progress_tick(&this->prog);
                    mpi_send_i(this->com, nex_slave, -1);
                    mpi_send_i(this->com, nex_slave, -1);
                }
                AssertMpi(MPI_Barrier(MPI_COMM_WORLD));
                //      progress_done(&this->prog);
            } else {
                /* non-master nodes */
                send_result();
                mpi_receive_i(this->com, 0, j);
                mpi_receive_i(this->com, 0, i);
                //      fprintf(stderr,"NODE %d receive move (%d %d)\n",this->com.getMyRank(),*j,*i);
                this->actual_row = (*i);
                this->actual_col = (*j);
                if((*j) < 0) {
                    AssertMpi(MPI_Barrier(MPI_COMM_WORLD));
                }
                return *j >= 0;
            }
            return 0;
#endif
        } else {
            /* no MPICC or only single node */
            (*j) = this->actual_col = newcol();
            (*i) = this->actual_row = this->count_row;
            /*    if (this->actual_col<0)
            progress_done(&this->prog);
            else {
              progress_tick(&this->prog);
              } */
            return this->actual_col >= 0;
        }
        return 0;
    }
};
}

#endif
