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
/*TODO: change to ! \file
 * \brief
 * Implements FcaTask2D.
 */
#include "gmxpre.h"

#include "task_2d.h"

#include "config.h"

#include <utility>

#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

#if GMX_LIB_MPI
int FcaTask2D::receive_result()
{
    int    slave;
    double data;
    int    i, j;
    slave = mpi_receive_request(this->com);
    mpi_receive_d(this->com, slave, &data);
    mpi_receive_i(this->com, slave, &i);
    mpi_receive_i(this->com, slave, &j);
    /* receives one bogus-result for cell (-1,-1) in beginning from every node, do not store... */
    if (i >= 0)
    {
        this->mat[i * this->M + j] = data;
        if (this->bSym)
        {
            this->mat[j * this->N + i] = data;
        }
    }
    return slave;
}

void FcaTask2D::send_result()
{
    mpi_send_request(this->com);
    mpi_send_d(this->com, 0, this->result);
    mpi_send_i(this->com, 0, this->actual_row);
    mpi_send_i(this->com, 0, this->actual_col);
}
#endif

FcaTask2D::FcaTask2D(const fca_utils::mpi &inCom, double* mat, const int M,
                     const int N, const bool bSym)
    : com(inCom)
{
    this->count_col  = 0;
    this->count_row  = 0;
    this->actual_row = -1;
    this->actual_col = -1;
    this->M          = M;
    this->N          = N;
    this->bSym       = bSym;
    this->mat        = mat;
    this->to_compute = nullptr;
    if (bSym && M != N)
    {
        fprintf(stdout, "init requesting a symmetric matrix with different dimensions %d x %d\n", M, N);
        exit(1);
        /* replaced this due to compatibility issues gmx3.3 vs gmx3.2
           gmx_fatal(FARGS,"init requesting a symmetric matrix with different dimensions %d x %d\n",M,N); */
    }
    if (com.getMyRank() == FCA_MASTER_RANK)
    {
        this->prog = FcaProgress("", (bSym ? N * (N + 1) / 2 : N * this->N),
                                 100,                   /* parallel tasks are very slow ... get more ticks to view */
                                 -com.getNumNodes() + 1 /* in the beginning there is a start tick for all nodes this is no real progress. */
                                 );
    }
    /* for sym: N*(N+1)/2 since diagonal elements are ticked, too */
    this->pairs      = nullptr;
    this->pair_count = -1;
    /* Initialize to values, even if they are not used */
    this->result            = -1.;
    this->compute_threshold = -1.;
}

FcaTask2D::FcaTask2D(const fca_utils::mpi &com, double* mat, const int M, const int N,
                     const bool bSym, const std::pair<int, int> pairs[], const int npairs)
    : FcaTask2D(com, mat, M, N, bSym)
{
    this->pairs      = pairs;
    this->pair_count = npairs;
}

void FcaTask2D::set_to_compute(real* to_compute, real threshold)
{
    this->to_compute        = to_compute;
    this->compute_threshold = threshold;
}

int FcaTask2D::finished()
{
    if (this->pairs)
    {
        return this->count_col >= this->pair_count;
    }
    else
    {
        return this->count_row >= this->M;
    }
}

int FcaTask2D::newcol()
{
    if (this->pairs)
    {
        if (!finished())
        {
            while (this->pairs[this->count_col].first <= 0)
            {
                this->count_col++;
                if (this->count_col >= this->pair_count)
                {
                    return -1;
                }
            }
            this->count_row = this->pairs[this->count_col].first;
            return this->pairs[(this->count_col)++].second;
        }
        else
        {
            return -1;
        }
    }
    else
    {
        if (((this->count_col > this->count_row) && this->bSym) || (this->count_col >= this->N))
        {
            this->count_col = 0;
            ++(this->count_row);
        }
        if (!finished())
        {
            if (this->to_compute)
            {
                if (this->to_compute[this->count_row + this->count_col * this->M] > this->compute_threshold)
                {
                    this->to_compute[this->count_row + this->count_col * this->M] = 0;
                    return (this->count_col)++;
                }
                else
                {
                    ++this->count_col;
                    return newcol();
                }
            }
            else //that is !to_compute
            {
                return (this->count_col)++;
            }
        }
        else
        {
            return -1;
        }
    }
}

void FcaTask2D::push_result(double res)
{
    this->result = res;
    if (this->com.getNumNodes() == 1)
    {
        /* no GMX_MPI or only single node */
        this->mat[this->actual_row * this->M + this->actual_col] = res;
        if (this->bSym)
        {
            this->mat[this->actual_col * this->N + this->actual_row] = res;
        }
    }
}

int FcaTask2D::next_index(int* i, int* j)
{
    if (this->com.getNumNodes() > 1)
    {
#if GMX_LIB_MPI
        if (this->com.getMyRank() == FCA_MASTER_RANK)
        {
            /*one job is for the master himself */
            int col = newcol();
            while (col >= 0)
            {
                int slave = receive_result();
                //          progress_tick(&this->prog); /* tick before the job is done */
                mpi_send_i(this->com, slave, col);
                mpi_send_i(this->com, slave, this->count_row);
                col = newcol();
            }
            for (int slave = 1; slave < this->com.getNumNodes(); slave++)
            {
                int nex_slave = receive_result();
                //	progress_tick(&this->prog);
                mpi_send_i(this->com, nex_slave, -1);
                mpi_send_i(this->com, nex_slave, -1);
            }
            fca_utils::AssertMpi(MPI_Barrier(MPI_COMM_WORLD), __FILE__, __LINE__);
            //      progress_done(&this->prog);
        }
        else
        {
            /* non-master nodes */
            send_result();
            mpi_receive_i(this->com, 0, j);
            mpi_receive_i(this->com, 0, i);
            //      fprintf(stdout,"NODE %d receive move (%d %d)\n",this->com.getMyRank(),*j,*i);
            this->actual_row = (*i);
            this->actual_col = (*j);
            if ((*j) < 0)
            {
                fca_utils::AssertMpi(MPI_Barrier(MPI_COMM_WORLD), __FILE__, __LINE__);
            }
            return *j >= 0;
        }
#endif
        return 0;
    }
    else
    {
        /* no GMX_MPI or only single node */
        (*j) = this->actual_col = newcol();
        (*i) = this->actual_row = this->count_row;
        return this->actual_col >= 0;
    }
}
