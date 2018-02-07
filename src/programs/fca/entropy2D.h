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
 * Declares FcaEntropy2D.
 */
#ifndef FCA_ENTROPY2D_H
#define FCA_ENTROPY2D_H

#include <memory>
#include <utility>

#include "gromacs/utility/real.h"

#include "entropy_planes.h"

class FcaEntropy2D
{
    int nbins;                         /* M bins : number of bins in bins and wbins */
    std::unique_ptr< real[] > bins;    /* MxM vector of count per bin */
    std::unique_ptr< real[] > p_quant; /* MxM */
    // not used ??? real *wbins;

    int ndata;                                         /* N data points */
    std::unique_ptr< std::pair<real, real>[] > eta;    /* 1xN data-point i is placed at position 0<=eta[i] <1 between bin binndx[i] and binndx[i]+1; */
    std::unique_ptr< std::pair<int, int>[] >   binndx; /* 1xN index of bin in which one finds data point i */

    /* a precomputed gauss kernel */
    int  ngauss;                     /* S points in function table for gauss */
    std::unique_ptr< real[] > gauss; /* 1xS function values for gauss */
    real dvol;                       /* */
    std::unique_ptr< real[] > p;     /* */

    void entropy2D_binning(const real sx[], const real sy[]);

    void entropy2D_convolution(const real kernel[], const real bins[], real result[]);

    real entropy2D_interpolate() const;

    void setup_bins(const real s[], real* minv, real* maxv, real* dx);

    public:
        /* S - entropy2D, gS - entropy2D gradient */

        FcaEntropy2D(const int ndata, const int nbins, const bool bAccurate);

        virtual ~FcaEntropy2D() { }

        real entropy2D(const real s1[], const real s2[]);

        real entropy2D_plane(const FcaPlane &plane) { return entropy2D(plane.getS1(), plane.getS2()); }
};

#endif
