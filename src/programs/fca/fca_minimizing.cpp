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
 * Implements FcaMaster.
 */
#include "gmxpre.h"

#include "fca_minimizing.h"

#include "config.h"

#include <cassert>
#include <cmath>

#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "task_2d.h"

void FcaMaster::rand_array_int(gmx::ThreeFry2x64<64> &rng,
                               const int n, int array[])
{
    for (int i = 0; i < n; i++)
    {
        array[i] = i;
    }

    int i = 0;
    while (i < 5 * n)
    {
        const int idx1 = random_int(rng, n);
        const int idx2 = random_int(rng, n);
        std::swap(array[idx1], array[idx2]);
        i += 1;
    }
}

void FcaMaster::update_icmatrix(const int nr, const std::pair<int, int>* moves)
{
    for (int i = 0, ct = 0; ct < nr; i++)
    {
        std::pair<int, int> localMove(0, 0);
        double              theta = 0;

        if (com.isMaster())
        {
            if (moves[i].first >= 0)
            {
                ct++;
                theta            = this->theta[moves[i].first + moves[i].second * dim];
                localMove.first  = moves[i].first;
                localMove.second = moves[i].second;
            }
            else
            {
                continue;
            }
        }
        else
        {
            //slave-nodes kriegen nur sinnvolle moves gesendet.
            ct += 1;
        }
#if GMX_LIB_MPI
        if (this->com.getNumNodes() > 1)
        {
            fca_utils::AssertMpi(MPI_Bcast(&localMove, 2, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
            fca_utils::AssertMpi(MPI_Bcast(&theta, 1, MPI_DOUBLE, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
        }
#endif
        const real costh      = cos(theta);
        const real sinth      = sin(theta);
        real     * ic1        = this->ic_matrix[localMove.first].get();
        real     * ic2        = this->ic_matrix[localMove.second].get();
        for (int j = 0; j < dim; j++)
        {
            real dum = ic1[j];
            ic1[j]   = costh * dum + sinth * ic2[j];
            ic2[j]   = costh * ic2[j] - sinth * dum;
        }

        if (com.isMaster())
        {
            const int a = localMove.first;
            const int b = localMove.second;
            theta = fmod(theta, fca_utils::PI_constant() / 2);
            if (theta > fca_utils::PI_constant() / 4)
            {
                theta = theta - fca_utils::PI_constant() / 2;
            }
            const real atheta = fabs(theta);
            for (int j = 0; j < dim; j++)
            {
                this->move[b + j * dim]       += atheta;
                this->move[j + a * dim]       += atheta;
                this->to_compute[a + j * dim] += atheta;
                this->to_compute[j + b * dim] += atheta;
            }
            this->move[b + a * dim] = 0;
        } //is Master
    }     // for nmoves
}         //update_fca

std::vector < std::pair < int, int>> FcaMaster::generate_moves(const int nr,
                                                               fca_utils::Log& /*logger*/)
{
    //computes the new MI-matrix and uses
    //moves*MI to find nr new moves
    //the algorithm works as follows. The 2xnr array moves and the 1xnr array val are filled with values and their indices from
    //MI*this->move.
    //min_val always points to the lowest entry in moves/val. a value from MI is compared with this entry.
    //If higher than the lowest entry, the lowest entry is replaced by the value from MI, and we search the 1xnr array val,
    //to find the now lowest entry, and set min_val to this entry.
    //Function tested with MATLAB.

    if (com.isMaster() == false)
    {
        return std::vector < std::pair < int, int>>();
    }

    // go through every entry of the MI-matrix
    for (int i = 1; i < dim; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (min_MI > this->MI[i + j * dim])
            {
                min_MI = this->MI[i + j * dim];
            }
        }
    }

    struct MoveState
    {
        std::pair<int, int> indexes;
        real                dum;
    };

    std::vector<MoveState> currentStates;
    currentStates.reserve(dim*(dim/2));

    for (int i = 1; i < dim; i++)
    {
        for (int j = 0; j < i; j++)
        {
            const real dum = (this->MI[i + j * dim] - min_MI) * this->move[j + i * dim];

            // if we really want to do this move...
            if (this->move[j + i * dim] > 0.2 && (this->MI[i + j * dim] + min_MI) > MI_MIN_CORR)
            {
                currentStates.emplace_back(MoveState {{i, j}, dum});
            }
        }
    }   //finished with MI-Matrix

    fca_utils::PSort::sort(currentStates.begin(), currentStates.end(), [](const MoveState &st1,
                                                                          const MoveState &st2)
                           {
                               return st1.dum > st2.dum;
                           });

    std::vector<bool> indexesUsed(dim, false);
    std::vector < std::pair < int, int>> nextMoves;
    nextMoves.reserve(dim);

    for (const MoveState &st : currentStates)
    {
        if (indexesUsed[st.indexes.first] == false
            && indexesUsed[st.indexes.second] == false)
        {
            nextMoves.emplace_back(st.indexes);
            indexesUsed[st.indexes.first]  = true;
            indexesUsed[st.indexes.second] = true;

            // we stop if all indexes are used in a pair
            if (static_cast<int>(nextMoves.size()) == dim/2)
            {
                break;
            }
            // we stop if nr pairs have been found
            if (static_cast<int>(nextMoves.size()) == nr)
            {
                break;
            }
        }
    }

    return nextMoves;
}

FcaMaster::FcaMaster(const int inDim, const int innframes,
                     std::vector< std::unique_ptr< real[] > > inProjx,
                     const fca_utils::mpi &inCom)
    : com(inCom),
      useLogInMinimization(false),
      verboseEnable(true),
      symmetricTask(true),
      plane(innframes),
      projx(std::move(inProjx)),
      MIsum(0)
{
    this->nframes = innframes;
    this->dim     = inDim;
    this->steps   = 0;
    this->bSync   = true;
    this->min_MI  = 100;
    // broadcast dim,nframes and projx
#if GMX_LIB_MPI
    if (com.getNumNodes() > 1)
    {
        {
            int nframesdim[2];
            if (com.isMaster())
            {
                nframesdim[0] = nframes;
                nframesdim[1] = dim;
            }
            fca_utils::AssertMpi(MPI_Bcast(nframesdim, 2, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
            if (com.isMaster() == false)
            {
                nframes = nframesdim[0];
                dim     = nframesdim[1];
            }
        }

        if (com.isMaster() == false)
        {
            assert(projx.size() == 0);
            projx.resize(nframes);
            for (int i = 0; i < nframes; i++)
            {
                projx[i].reset(new real[dim]);
            }
        }

        const MPI_Datatype realtype = (std::is_same<real, float>::value ? MPI_FLOAT : MPI_DOUBLE);
        for (int i = 0; i < nframes; i++)
        {
            fca_utils::AssertMpi(MPI_Bcast(projx[i].get(), dim, realtype, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
        }
    }
#endif

    if (com.isMaster())
    {
        this->to_compute.reset(new real[dim*dim]());
        this->move = this->to_compute.get();
        this->theta.reset(new double[dim*dim]());
        this->MI.reset(new double[dim*dim]());
        for (int i = 0; i < dim; i++)
        {
            for (int k = 0; k < dim; k++)
            {
                this->to_compute[i + k * dim] = 1.0; /*setze move und to_compute auf 1 */
            }
        }
    }

    /*  allocate ic_matrix */
    this->ic_matrix.reset(new std::unique_ptr< real[] >[dim]);
    this->icx.reset(new std::unique_ptr< real[] >[dim]);
    this->S1.reset(new real[dim]);
    for (int i = 0; i < dim; i++)
    {
        this->icx[i].reset(new real[nframes]);
        this->ic_matrix[i].reset(new real[dim]);
        for (int k = 0; k < dim; k++)
        {
            this->ic_matrix[i][k] = 0;
        }
        this->ic_matrix[i][i] = 1;
    }
}

void FcaMaster::dump_ic_matrix(FILE* out) const
{
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            fprintf(out, "%e ", this->ic_matrix[i][j]);
        }
        fprintf(out, "\n");
    }
}

void FcaMaster::read_ic_matrix(const char inFilename[])
{
    FILE* in = gmx_ffopen(inFilename, "r");
    assert(in);

    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            int nread = fca_utils::fscanf_floating_point(in, &this->ic_matrix[i][j]);
            if (nread < 1)
            {
                fprintf(stdout, "FATAL: file does not contain enough elements to read (%dx%d)-matrix\n", dim, dim);
                exit(1);
            }
        }
    }
    fclose(in);

    this->bSync = false;
}

void FcaMaster::broadcast_ic_matrix()
{
#if GMX_LIB_MPI
    fca_utils::AssertMpi(MPI_Bcast(&this->bSync, sizeof(bool), MPI_BYTE, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    if (!this->bSync && (this->com.getNumNodes() > 1))
    {
        for (int i = 0; i < dim; i++)
        {
            fca_utils::AssertMpi(MPI_Bcast(this->ic_matrix[i].get(), sizeof(real) * dim, MPI_BYTE, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
        }
    }
#endif
    this->bSync = true;
}

/*TODO: change to ! \brief
 * Returns zero if no move has been possible
 */
int FcaMaster::minimization(const int nr_moves, fca_utils::Log &logger)
{
    compute_MI_matrix(this->to_compute.get(), nullptr);
    this->steps += 1;

    std::vector < std::pair < int, int>> next_moves;

    if (com.isMaster())
    {
        if (logger.MI)
        {
            dump_MIsum(logger.MI);
        }
        if (useLogInMinimization && this->steps % LOGFRQ == 0)
        {
            if (logger.move_matrix)
            {
                fca_utils::dump_sqrmatrixr(logger.move_matrix, this->to_compute.get(), this->dim);
            }
            if (logger.ic_matrix)
            {
                fprintf(logger.ic_matrix, "================= step %d ============================\n", this->steps);
                dump_ic_matrix(logger.ic_matrix);
            }
            if (logger.MI_matrix)
            {
                fprintf(logger.MI_matrix, "================= step %d ============================\n", this->steps);
                write_MI_matrix(logger.MI_matrix);
            }
        }
        next_moves = generate_moves(nr_moves, logger);
        assert(int(next_moves.size()) <= nr_moves);
    }

    int actual_moves = static_cast<int>(next_moves.size());
#if GMX_LIB_MPI
    fca_utils::AssertMpi(MPI_Bcast(&actual_moves, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
#endif
    if (actual_moves > 0)
    {
        FcaTask2D task2D(this->com, this->theta.get(), this->dim, this->dim, symmetricTask, next_moves.data(), int(next_moves.size()));
        if (useLogInMinimization && com.isMaster())
        {
            fprintf(stdout, "NODE %d, starting plane rotations...\n", this->com.getMyRank());
        }
        // auch in fca-struct rein ?
        FcaEntropy entr(this->nframes, 200);

        int        p1, p2;
        while (task2D.next_index(&p1, &p2))
        {
            if (p1 != p2)
            {
                real theta;
                this->plane.init_plane(this->projx.data(), this->ic_matrix[p1].get(), this->ic_matrix[p2].get(), this->dim);
                theta = this->plane.minimize_plane(&entr);
                task2D.push_result(theta);
            }
            else
            {
                fprintf(stdout, "NODE %d: warum soll ich plane (%d,%d) rotieren ?? \n", this->com.getMyRank(), p1, p2);
                task2D.push_result(0);
            }
        }
        update_icmatrix(actual_moves, next_moves.data());
    } // if (actual moves >0)
    return actual_moves > 0;
}

void FcaMaster::project_ic_signal()
{
    if (verboseEnable)
    {
        fprintf(stdout, "project signal to ics\n");
    }
#if GMX_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int i = 0; i < dim; i++)
    {
        for (int k = 0; k < nframes; k++)
        {
            icx[i][k] = 0;
            for (int j = 0; j < dim; j++)
            {
                icx[i][k] += ic_matrix[i][j] * projx[k][j];
            }
        }
    }
}

void FcaMaster::compute_S1()
{
    // computes S1 for all ic-modes
    // non-parallel version
    // ACHTUNG, wenn parallelisiert darf es nicht mehr von
    // write_fca_valfile aufgerufen werden
    // real *to_compute=this->to_compute;
    // t_mpi *com=this->com;
#if GMX_OPENMP
    std::unique_ptr<std::unique_ptr<FcaEntropy>[]> entr1DforThreads(new std::unique_ptr<FcaEntropy>[omp_get_max_threads()]);
#else
    FcaEntropy entr1D(nframes, 200);
#endif

    this->MIsum = 0;

    if (verboseEnable)
    {
        fprintf(stdout, "single S 0-%d  ", dim - 1);
        fprintf(stdout, "berechne 1D entropies...\n");
    }

    real localMIsum = 0;
#if GMX_OPENMP
#pragma omp parallel for schedule(static) default(shared) reduction(+:localMIsum)
#endif
    for (int i = 0; i < dim; i++)
    {
#if GMX_OPENMP
        if (entr1DforThreads[omp_get_thread_num()] == nullptr)
        {
            entr1DforThreads[omp_get_thread_num()].reset(new FcaEntropy(nframes, 200));
        }

        entr1DforThreads[omp_get_thread_num()]->entropy_1D(&S1[i], nullptr, icx[i].get(), nullptr);
#else
        entr1D.entropy_1D(&S1[i], nullptr, icx[i].get(), nullptr);
#endif
        if (verboseEnable)
        {
            fprintf(stdout, "%e ", S1[i]);
        }

        localMIsum += S1[i];
    }

    this->MIsum += localMIsum;
}

void FcaMaster::compute_MI_matrix(real* to_compute, gmx::ThreeFry2x64<64>* rng)
{
    std::unique_ptr< int[] > permute(new int[this->nframes]());
    std::unique_ptr< std::unique_ptr< real[] >[] > perm_data(new std::unique_ptr< real[] >[this->dim]());
    if (rng)
    {
        for (int d = 0; d < this->dim; d++)
        {
            perm_data[d].reset(new real[this->nframes]());
        }
    }
    const real computeAgain = 0.2;
    project_ic_signal();
    if (rng)
    {
        for (int i = 0; i < dim; i++)
        {
            rand_array_int(*rng, nframes, permute.get());
            for (int j = 0; j < nframes; j++)
            {
                perm_data[i][j] = icx[i][permute[j]];
            }
        }
    }
    else
    {
        compute_S1();
    }

    if (verboseEnable && com.isMaster())
    {
        fprintf(stdout, "\n");
        fprintf(stdout, "berechne 2D entropies...");
    }

    FcaTask2D    task2D(com, MI.get(), dim, dim, symmetricTask);
    task2D.set_to_compute(to_compute, computeAgain);
    FcaEntropy2D entr2D(nframes, 100, rng != nullptr);

    int          p1, p2;
    while (task2D.next_index(&p1, &p2))
    {
        if (p1 != p2)
        {
            if (rng)
            {
                task2D.push_result(entr2D.entropy2D(perm_data[p1].get(), perm_data[p2].get()) - entr2D.entropy2D(icx[p1].get(), icx[p2].get()));
            }
            else
            {
                task2D.push_result(S1[p1] + S1[p2] - entr2D.entropy2D(icx[p1].get(), icx[p2].get()));
            }
        }
        else
        {
            task2D.push_result(0);
        }
    }
}

void FcaMaster::sort_modes(real* amplitude, real* anharmonicity,
                           const real WEIGHT_AMP, const real WEIGHT_ANH)
{
    struct qsort_fca
    {
        real                      amp;
        real                      anh;
        std::unique_ptr< real[] > ic;
    };

    std::unique_ptr< qsort_fca[] > qsort_data(new qsort_fca[dim]);
    for (int i = 0; i < dim; i++)
    {
        qsort_data[i].amp = amplitude[i];
        qsort_data[i].anh = anharmonicity[i];
        qsort_data[i].ic  = std::move(ic_matrix[i]);
    }
    fca_utils::PSort::sort(&qsort_data[0], &qsort_data[dim], [&](const qsort_fca &a, const qsort_fca &b) {
                               const gmx_bool res = (WEIGHT_AMP * b.amp + WEIGHT_ANH * log(1e-9 + std::abs(b.anh)))
                                   > ((a.amp) * WEIGHT_AMP + log(1e-9 + std::abs(a.anh)) * WEIGHT_ANH);
                               return (res);
                           });

    for (int i = 0; i < dim; i++)
    {
        amplitude[i]     = qsort_data[i].amp;
        anharmonicity[i] = qsort_data[i].anh;
        ic_matrix[i]     = std::move(qsort_data[i].ic);
    }
    bSync = FALSE;
}

void FcaMaster::compute_mode_amplitude(real* amplitude) const
{
    const real invnframes = 1.0 / nframes;
#if GMX_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int k = 0; k < dim; k++)
    {
        real amplitude_k = 0;
        for (int j = 0; j < nframes; j++)
        {
            amplitude_k += fca_utils::squareof(icx[k][j]);
        }
        amplitude_k *= invnframes;
        amplitude[k] = amplitude_k;
    }
}

void FcaMaster::compute_mode_anharmonicity(const real* amplitude,
                                           real      * anharm) const
{
    const real logsp = log(2 * fca_utils::PI_constant());
    for (int i = 0; i < dim; i++)
    {
        anharm[i] = 0.5 * (1 + logsp + log(amplitude[i])) - S1[i];
    }
}
