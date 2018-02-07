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
 * Declares FcaMaster.
 */
#ifndef FCA_MINIMIZING_H
#define FCA_MINIMIZING_H

#include <memory>
#include <utility>
#include <vector>

#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/utility/real.h"

#include "entropy2D.h"
#include "utils.h"

namespace gmx
{

class FcaMaster
{
    static constexpr real     MI_MIN_CORR = -0.5;
    static constexpr int      LOGFRQ      = 50;

    int                       dim;
    int                       nframes;

    bool                      bSync; /* TRUE if ic_matrix is the same on all nodes */

    /* References */
    const fca_utils::mpi     &com;
    const bool                useLogInMinimization;
    const bool                verboseEnable;
    const bool                symmetricTask;

    /* Owned objects: */
    /* ...used on all nodes */
    std::unique_ptr< std::unique_ptr< real[] >[] > ic_matrix;
    FcaPlane plane;
    std::vector< std::unique_ptr< real[] > >       projx;
    std::unique_ptr< std::unique_ptr< real[] >[] > icx;
    std::unique_ptr< real[] > S1;
    real MIsum;
    /* ...only used on MASTER */
    std::unique_ptr< double[] > theta;
    // double, damit wir beim MPI-stuff keine probleme bekommen
    std::unique_ptr< double[] > MI;
    //  move und to_compute zeigen auf den selben speicherbereich!
    std::unique_ptr< real[] >   to_compute;
    // i+k*dim  k<i --> to_compute / k>i --> move
    real * move;
    int    steps;
    double min_MI; /* MI shifts systematically with nframes -- to correct for this use min_MI*/

    static int random_int(ThreeFry2x64<64> &rng, int maxv) { return UniformIntDistribution<int>(0, maxv-1) (rng); }

    /*! \brief
     * Any permutaion of N numbers can be reached by pairwise permutations
     * 5*N seems sufficient to randomize these indices
     */
    static void rand_array_int(ThreeFry2x64<64> &rng, const int n, int array[]);

    /*! \brief
     * Returns zero if no move has been possible
     */
    void update_icmatrix(const int nr, const std::pair<int, int>* moves);

    std::vector < std::pair < int, int>> generate_moves(const int nr, fca_utils::Log& /*logger*/);

    void dump_MIsum(FILE* out) const { fprintf(out, "%e\n", this->MIsum); }

    public:
        int getNFrames() const { return nframes; }
        std::unique_ptr< real[] >* getIcx() const { return icx.get(); }
        int getDim() const { return dim; }
        real* getS1() { return S1.get(); }
        double* getMI() { return MI.get(); }
        real getMIsum() const { return MIsum; }

        std::unique_ptr< real[] >* getIcMatrix() { return ic_matrix.get(); }

        FcaMaster(const int inDim, const int innframes, std::vector< std::unique_ptr< real[] > > inProjx, const fca_utils::mpi &inCom);

        ~FcaMaster() { }

        void write_MI_matrix(FILE* out) { fca_utils::dump_sqrmatrix(out, this->MI.get(), this->dim); }

        void dump_ic_matrix(FILE* out) const;

        void read_ic_matrix(const char inFilename[]);

        void broadcast_ic_matrix();

        int minimization(const int nr_moves, fca_utils::Log &logger);

        void project_ic_signal();

        /*! \brief
         * Input: fca - object algo - Shwartz or basic
         */
        void compute_S1();

        /*! \brief
         * Compute a matrix of pairwise correlations:
         * in theory: I(i,j) = S(i) + S(j) - S( i,j)  EQ(1)
         * however, 1D and 2D histograms have different effect on data ---> spurious correlations especially if I(i,j) approx 0
         * --> instead compute I( i, j) = S(i,j) - S*(i,j) where S* is computed from  data-series X(k),Y(k) by computing
         *                       S*(X,Y)= S( X(k), Y(P(k)) ), where P(k) is a random permutation of the indices k
         *                       INPUT tFCA -- fca object
         *                       to_compute -- which dimensions to compute
         *                       rng   -- if nullptr EQ(1) is used, otherwise rng is used to get permuation P(k) of indices
         */
        void compute_MI_matrix(real* to_compute, ThreeFry2x64<64>* rng);

        void sort_modes(real* amplitude, real* anharmonicity,
                        const real WEIGHT_AMP, const real WEIGHT_ANH);

        /*! \brief
         * Computes the amplitude of projections on ic-modes: icx
         * RETURN:   amplitude[0..dim-1]
         */
        void compute_mode_amplitude(real* amplitude) const;

        /*! \brief
         * Computes the negentropy J[x_i] =  Hgauss[x_i]-H[x_i]
         * INPUT
         *  S1  entropy of mode   [0..dim-1]
         *  amplitude .. of mode  [0..dim-1]
         * RETURN
         *  anharm[0..dim-1]
         */
        void compute_mode_anharmonicity(const real* amplitude, real* anharm) const;
};

} //gmx namespace

#endif
