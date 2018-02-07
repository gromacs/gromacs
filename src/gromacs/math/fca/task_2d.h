/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares FcaProgress and FcaTask2D.
 */
#ifndef FCA_TASK_2D_H
#define FCA_TASK_2D_H

#include "config.h"

#include <cstdio>

#include <string>
#include <utility>

#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

#include "utils.h"

namespace gmx
{

class FcaProgress
{
    int         pos;
    int         maxticks;
    int         shown_ticks;
    int         full;
    std::string msg;
    FILE      * out;

    void progress_show()
    {
        const int now  = lrint(1.0 * this->pos / this->maxticks * full);

        if (now > this->shown_ticks)
        {
            this->shown_ticks = now;
            progress_resetline();
            fprintf(this->out, "%s", this->msg.c_str());
            fprintf(this->out, " ");
            for (int i = 0; i < now; i++)
            {
                fprintf(this->out, "*");
            }
            for (int i = now; i < full; i++)
            {
                fprintf(this->out, ".");
            }
        }
    }

    void progress_resetline() { fprintf(this->out, "\r"); }

    void progress_tick() { this->pos += 1; progress_show(); }

    void progress_forward(const int plus) { this->pos += plus; progress_show(); }

    void progress_done() { fprintf(this->out, "\n"); }

    public:
        FcaProgress(FILE* inOut, const std::string &inMsg, const int inMaxticks,
                    const int inFull, const int inPos)
            : pos(inPos),
              maxticks(inMaxticks),
              shown_ticks(-1),
              full(inFull),
              msg(inMsg),
              out(inOut) { }

        FcaProgress(const std::string &inMsg, const int inMaxticks,
                    const int inFull, const int inPos)
            : FcaProgress(stdout, inMsg, inMaxticks, inFull, inPos) { }

        explicit FcaProgress() : FcaProgress(stdout, "", 0, 0, 0) { }

        FcaProgress(const FcaProgress &)            = default;
        FcaProgress &operator=(const FcaProgress &) = default;

        virtual ~FcaProgress() { }
};

class FcaTask2D
{
    static const int               __TAG       = 42;
    static const int               REQUEST_TAG = 44;

    const fca_utils::mpi          &com;
    double                       * mat;
    int                            N, M;
    int                            count_row;
    int                            count_col;
    int                            actual_row;
    int                            actual_col;
    int                            bSym;
    double                         result;
    const std::pair<int, int>    * pairs;      /* if set, only compute these pairs */
    int                            pair_count;
    real                         * to_compute; /* compute value only if entry in to_compute larger than "compute_threshold" */
    real                           compute_threshold;
    FcaProgress                    prog;

#if GMX_LIB_MPI
    static void _receive(void* data, int nr, MPI_Datatype type, int address, int tag, MPI_Comm comm, MPI_Status* stat) { fca_utils::AssertMpi(MPI_Recv(data, nr, type, address, tag, comm, stat), __FILE__, __LINE__); }

    static void _send(void* dat, int nr, MPI_Datatype type, int address, int tag, MPI_Comm comm) { fca_utils::AssertMpi(MPI_Send(dat, nr, type, address, tag, comm), __FILE__, __LINE__); }

    static void mpi_send_it(const fca_utils::mpi & /*mpiInfo*/, int address, int dat, int tag) { _send(&dat, 1, MPI_INT, address, tag, MPI_COMM_WORLD); }

    static void mpi_send_i(const fca_utils::mpi &mpiInfo, int address, int dat) { mpi_send_it(mpiInfo, address, dat, __TAG); }

    static void mpi_send_d(const fca_utils::mpi & /*mpiInfo*/, int address, double dat) { _send(&dat, 1, MPI_DOUBLE, address, __TAG, MPI_COMM_WORLD); }

    static void mpi_receive_d(const fca_utils::mpi & /*mpiInfo*/, int address, double* data) { MPI_Status stat; _receive(data, 1, MPI_DOUBLE, address, __TAG, MPI_COMM_WORLD, &stat); }

    static void mpi_receive_i(const fca_utils::mpi & /*mpiInfo*/, int address, int* data) { MPI_Status stat; _receive(data, 1, MPI_INT, address, __TAG, MPI_COMM_WORLD, &stat); }

    static int mpi_receive_request(const fca_utils::mpi & /*mpiInfo*/)
    {
        int        data;
        MPI_Status stat;
        _receive(&data, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &stat);
        return data;
    }

    static void mpi_send_request(const fca_utils::mpi &mpiInfo)
    {
        mpi_send_it(mpiInfo, 0, mpiInfo.getMyRank(), REQUEST_TAG);
    }
#endif

    int finished();

    int newcol();

    public:

#if GMX_LIB_MPI
        int receive_result();

        void send_result();
#endif

        FcaTask2D(const fca_utils::mpi &inCom, double* mat, const int M, const int N, const bool bSym);

        FcaTask2D(const fca_utils::mpi &com, double* mat, const int M, const int N, const bool bSym, const std::pair<int, int> pairs[], const int npairs);

        virtual ~FcaTask2D() { }

        FcaTask2D(const FcaTask2D &)            = delete;
        FcaTask2D &operator=(const FcaTask2D &) = delete;

        void set_to_compute(real* to_compute, real threshold);

        void push_result(double res);

        int next_index(int* i, int* j);
};

} //gmx namespace

#endif
