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
#ifndef FCA_ENTROPY2D_H
#define FCA_ENTROPY2D_H

#include <cmath>

#include <memory>

#include "entropy_planes.h"
#include "task_2d.hpp"
#include "utils.hpp"
//#include "fca_minimizing.h" //for the tAlgoSelect

namespace FCA
{

class Entropy2D
{
    int nbins;                         /*M bins : number of bins in bins and wbins */
    std::unique_ptr< real[] > bins;    /*MxM vector of count per bin */
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

    void entropy2D_binning(const real sx[], const real sy[])
    {
        /*  the data s is binned (number of bins: this->nbins )
           the result is stored in entr
           eta -- interpolation distance to next left bin-edge 0<eta<1  eta=s[i]-binndx[i]
           binndx[i] -- in which bin is s[i]
           bins[j] -- the result of binning.
         */
        real xmax = FCA::utils::limits_min<real>();
        real xmin = FCA::utils::limits_max<real>();
        real ymax = FCA::utils::limits_min<real>();
        real ymin = FCA::utils::limits_max<real>();
        real dx, dy;

        setup_bins(sx, &xmin, &xmax, &dx);
        setup_bins(sy, &ymin, &ymax, &dy);
        this->dvol = dx * dy;
        /* initialize bins with zero */
        const int sqrenbins = utils::squareof(this->nbins);
        for (int i = 0; i < sqrenbins; i++)
        {
            this->bins[i] = 0.0;
        }

        /* bin the data */
        for (int i = 0; i < this->ndata; i++)
        {
            const int  xhash    = lrint(floor((sx[i] - xmin) / dx));
            const int  yhash    = lrint(floor((sy[i] - ymin) / dy));
            const real xeta     = (sx[i] - xmin) / dx - xhash;
            const real yeta     = (sy[i] - ymin) / dy - yhash;
            this->binndx[i].first                              = xhash;
            this->binndx[i].second                             = yhash;
            this->eta[i].first                                 = 1.0 - xeta;
            this->eta[i].second                                = 1.0 - yeta;
            this->bins[xhash + yhash * this->nbins]           += (1.0 - xeta) * (1.0 - yeta);
            this->bins[xhash + 1 + yhash * this->nbins]       += xeta * (1.0 - yeta);
            this->bins[xhash + (yhash + 1) * this->nbins]     += (1.0 - xeta) * yeta;
            this->bins[xhash + 1 + (1 + yhash) * this->nbins] += xeta * yeta;
        }
    }

    void entropy2D_convolution(const real kernel[], const real bins[], real result[])
    {
        /* this computes a convolution of bins with kernel
           could be replaced by fft convolution -- will that be faster ?
           now: 200 bins, ca 35 ngauss 200*35 operations: 7000 (direkt) vs 1035 (fft) but more overhead */
        const int  m     = floor(this->ngauss / 2.0);
        const real invn  = 1.0 / this->ndata;
        for (int x = 0; x < nbins; x++)
        {
            for (int y = 0; y < nbins; y++)
            {
                result[x + y * nbins] = 0;
                const int xgstart     = this->ngauss - (nbins - x) - m + 1; //bug-fix thanks to Shun Sakuraba 21/8/08
                const int xgstopp     = x + 1 + m;
                const int ygstart     = this->ngauss - (nbins - y) - m + 1; //bug-fix thanks to Shun Sakuraba 21/8/08
                const int ygstopp     = y + 1 + m;
                for (int i = 0 > xgstart ? 0 : xgstart; i < (this->ngauss < xgstopp ? this->ngauss : xgstopp); i++)
                {
                    for (int k = 0 > ygstart ? 0 : ygstart; k < (this->ngauss < ygstopp ? this->ngauss : ygstopp); k++)
                    {
                        result[x + y * nbins] += bins[x + m - i + (y + m - k) * nbins] * kernel[k] * kernel[i] * invn;
                    }
                }
            }
        }
    }

    real entropy2D_interpolate() const
    {
        real       sum           = 0.0;
        const real invdvol       = 1.0 / this->dvol;

        const int  nbins = this->nbins;
        for (int i = 0; i < this->ndata; i++)
        {
            const int  xhash = this->binndx[i].first;
            const int  yhash = this->binndx[i].second;
            const real xa    = this->eta[i].first;
            const real xb    = 1 - xa;
            const real ya    = this->eta[i].second;
            const real yb    = 1 - ya;
            const real p     = (xa * ya * this->p_quant[xhash + yhash * nbins] +
                                xb * ya * this->p_quant[xhash + 1 + yhash * nbins] +
                                xa * yb * this->p_quant[xhash + (yhash + 1) * nbins] +
                                xb * yb * this->p_quant[xhash + 1 + (yhash + 1) * nbins]) *
                invdvol;
            sum += -logf(p);
        }
        return sum * 1.0 / this->ndata;
    }

    public:
        /* S - entropy2D, gS - entropy2D gradient */

        Entropy2D(const int ndata, const int nbins, const bool bAccurate)
        {
            this->nbins = nbins;
            this->ndata = ndata;

            this->bins.reset(new real[nbins * nbins]());
            this->p_quant.reset(new real[nbins * nbins]());
            this->eta.reset(new std::pair<real, real>[ndata]());
            this->binndx.reset(new std::pair<int, int>[ndata]());

            real sigma;
            if (bAccurate)                                               // 1D gauss kernel of width sig
            {
                sigma        = 1 * (1 + 2 * log(ceil(20000.0 / ndata))); /* fuer 2000 Punkte ist 5 zu niedrig, aber bei 20000 ist es super */
                this->ngauss = ceil(sigma * 14);                         /* mit *14 ist es genauer, aber so reicht es eigentlich. cumsum(gradS) weicht zwar ein bisschen von S ab, aber hat selben hoch+tief punkte. */
            }
            else
            {
                sigma        = 1 * (1 + log10(ceil(20000.0 / ndata))); /* fuer 2000 Punkte ist 5 zu niedrig, aber bei 20000 ist es super */
                this->ngauss = ceil(sigma * 4);                        /* mit *14 ist es genauer, aber so reicht es eigentlich. cumsum(gradS) weicht zwar ein bisschen von S ab, aber hat selben hoch+tief punkte. */
            }

            this->gauss.reset(new real[this->ngauss]());

            for (int i = 0; i < this->ngauss; i++)
            {
                this->gauss[i] = 1.0 / sqrt(2.0 * M_PI) / sigma * exp(-utils::squareof(1.0 / sigma * (i - floor(this->ngauss / 2.0))) / 2);
            }
        }

        virtual ~Entropy2D()
        {
        }


        void setup_bins(const real s[], real* minv, real* maxv, real* dx)
        {
            /* find maxima and minima of data */
            for (int i = 0; i < this->ndata; i++)
            {
                (*maxv) = s[i] > (*maxv) ? s[i] : (*maxv);
                (*minv) = s[i] < (*minv) ? s[i] : (*minv);
            }

            /* compute bin-width dx */
            (*dx) = ((*maxv) - (*minv)) / (1.0 * this->nbins);

            /* correct for extra-space at the ends due to gauss-smoothing */
            /* correct to accomodate tails of gaussian-kernel */
            (*maxv) += ceil(this->ngauss / 2.0) * (*dx);
            (*minv) -= ceil(this->ngauss / 2.0) * (*dx);
            (*dx)    = ((*maxv) - (*minv)) / (1.0 * this->nbins);
        }

        real entropy2D(const real s1[], const real s2[])
        {
            entropy2D_binning(s1, s2);
            entropy2D_convolution(this->gauss.get(), this->bins.get(), this->p_quant.get());
            return entropy2D_interpolate();
        }

        real entropy2D_plane(const Plane &plane)
        {
            return entropy2D(plane.getS1(), plane.getS2());
        }
};
}

#endif
