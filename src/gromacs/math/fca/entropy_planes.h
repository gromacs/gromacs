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
 * Declares FcaEntropy and FcaPlane.
 */
#ifndef FCA_ENTROPY_PLANES_H
#define FCA_ENTROPY_PLANES_H

#include <cmath>

#include <memory>

#include "gromacs/math/vec.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class FcaEntropy
{
    int nbins;
    std::unique_ptr< real[] > bins;
    std::unique_ptr< real[] > wbins;
    std::unique_ptr< real[] > eta; /* data-point i is placed at position 0<=eta[i] <1 between bin binndx[i] and binndx[i]+1; */
    std::unique_ptr< int[] >  binndx;
    int ndata;
    std::unique_ptr< real[] > gauss;
    std::unique_ptr< real[] > gaussprime;
    int  ngauss;
    real dx;
    std::unique_ptr< real[] > p_quant;
    std::unique_ptr< real[] > phi_quant;
    std::unique_ptr< real[] > p;
    std::unique_ptr< real[] > phi;
    std::unique_ptr< real[] > Fquant;
    std::unique_ptr< real[] > Hgradl;

    void entropy_binning(const real s[]);

    void entropy_convolution(const real kernel[], const real bins[], real result[]) const;

    real entropy_interpolate_p();

    void entropy_interpolate_phi();

    real entropy_grad(const real gfun[]);

    public:
        FcaEntropy(const int inNdata, const int inNbins);

        virtual ~FcaEntropy() { }

        /*! \brief
         * results:
         * S -- resulting 1D entropy
         * gS -- nullptr or gradient of entropy
         * input :
         * s -- data-series of length this->ndata
         * g -- nullptr or ...
         * entr -- entropy object
         */
        real entropy_1D(real* S, real* gS, const real s[], const real g[]);
};

class FcaPlane
{
    std::unique_ptr< real[] > proj1; /* projections for theta=0; */
    std::unique_ptr< real[] > proj2;
    std::unique_ptr< real[] > s1;    /* projections for rotation in plane by theta */
    std::unique_ptr< real[] > s2;    /* s1,s2 for the normal projection */
    std::unique_ptr< real[] > g1;    /* for gradient */
    std::unique_ptr< real[] > g2;
    real theta;
    real S;
    real gS;
    int  ndata;

    public:
        explicit FcaPlane(const int inNdata);

        virtual ~FcaPlane() { }

        const real* getS1() const { return s1.get(); }

        const real* getS2() const { return s2.get(); }

        void init_plane(const std::unique_ptr< real[] > projx[], const real a1[], const real a2[], const int dim);

        void rotate_plane(const real theta);

        void entropy_plane(real* S, real* gS, FcaEntropy* entr);

        real minimize_plane(FcaEntropy* entr);

        static real SIGN(const real a, const real b)
        {
            return ((b) > 0.0 ? std::abs(a) : -std::abs(a));
        }

        static void MOV3(real &a, real &b, real &c, const real d, const real e, const real f)
        {
            (a) = (d);
            (b) = (e);
            (c) = (f);
        }

        real dbrent_plane(const real ax, const real bx, const real cx, FcaEntropy* entr, const real tol, real* xmin);
};

} //gmx namespace

#endif
