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
#include "gmxpre.h"

#include "entropy2D.h"

#include <cmath>

#include <memory>
#include <utility>

#include "gromacs/utility/real.h"

#include "utils.h"
//#include "fca_minimizing.h" //for the tAlgoSelect

namespace gmx
{

void FcaEntropy2D::entropy2D_binning(const real sx[], const real sy[])
{
    /*  the data s is binned (number of bins: this->nbins )
        the result is stored in entr
        eta -- interpolation distance to next left bin-edge 0<eta<1  eta=s[i]-binndx[i]
        binndx[i] -- in which bin is s[i]
        bins[j] -- the result of binning.
     */
    real xmax = fca_utils::limits_min<real>();
    real xmin = fca_utils::limits_max<real>();
    real ymax = fca_utils::limits_min<real>();
    real ymin = fca_utils::limits_max<real>();
    real dx, dy;

    setup_bins(sx, &xmin, &xmax, &dx);
    setup_bins(sy, &ymin, &ymax, &dy);
    this->dvol = dx * dy;
    /* initialize bins with zero */
    const int sqrenbins = fca_utils::squareof(this->nbins);
    for (int i = 0; i < sqrenbins; i++)
    {
        this->bins[i] = 0.0;
    }

    /* bin the data */
    for (int i = 0; i < this->ndata; i++)
    {
        const int  xhash    = lrint(std::floor((sx[i] - xmin) / dx));
        const int  yhash    = lrint(std::floor((sy[i] - ymin) / dy));
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

void FcaEntropy2D::entropy2D_convolution(const real kernel[], const real bins[], real result[])
{
    /* this computes a convolution of bins with kernel
       could be replaced by fft convolution -- will that be faster ?
       now: 200 bins, ca 35 ngauss 200*35 operations: 7000 (direkt) vs 1035 (fft) but more overhead */
    const int  m     = std::floor(this->ngauss / 2.0);
    const real invn  = 1.0 / this->ndata;
    for (int y = 0; y < nbins; y++)
    {
        for (int x = 0; x < nbins; x++)
        {
            result[x + y * nbins] = 0;
            const int xgstart     = std::max(0, this->ngauss - (nbins - x) - m + 1); //bug-fix thanks to Shun Sakuraba 21/8/08
            const int xgstopp     = std::min(x + 1 + m, this->ngauss);
            const int ygstart     = std::max(0, this->ngauss - (nbins - y) - m + 1); //bug-fix thanks to Shun Sakuraba 21/8/08
            const int ygstopp     = std::min(y + 1 + m, this->ngauss);
            for (int k = ygstart; k < ygstopp; k++)
            {
                real        sum_xy = 0;
                const real* pbins  = &bins[x + m + (y + m - k) * nbins];
                for (int i = xgstart; i < xgstopp; i++)
                {
                    // Apply  * kernel[i] * invn after the loop
                    sum_xy += pbins[-i] * kernel[i];
                }
                result[x + y * nbins] += sum_xy * kernel[k] * invn;
            }
        }
    }
}

real FcaEntropy2D::entropy2D_interpolate() const
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

FcaEntropy2D::FcaEntropy2D(const int ndata, const int nbins, const bool bAccurate)
{
    this->nbins = nbins;
    this->ndata = ndata;

    this->bins.reset(new real[nbins * nbins]());
    this->p_quant.reset(new real[nbins * nbins]());
    this->eta.reset(new std::pair<real, real>[ndata]());
    this->binndx.reset(new std::pair<int, int>[ndata]());

    real sigma;
    if (bAccurate)                                                         // 1D gauss kernel of width sig
    {
        sigma        = 1 * (1 + 2 * std::log(std::ceil(20000.0 / ndata))); /* fuer 2000 Punkte ist 5 zu niedrig, aber bei 20000 ist es super */
        this->ngauss = std::ceil(sigma * 14);                              /* mit *14 ist es genauer, aber so reicht es eigentlich. cumsum(gradS) weicht zwar ein bisschen von S ab, aber hat selben hoch+tief punkte. */
    }
    else
    {
        sigma        = 1 * (1 + log10(std::ceil(20000.0 / ndata))); /* fuer 2000 Punkte ist 5 zu niedrig, aber bei 20000 ist es super */
        this->ngauss = std::ceil(sigma * 4);                        /* mit *14 ist es genauer, aber so reicht es eigentlich. cumsum(gradS) weicht zwar ein bisschen von S ab, aber hat selben hoch+tief punkte. */
    }

    this->gauss.reset(new real[this->ngauss]());

    for (int i = 0; i < this->ngauss; i++)
    {
        this->gauss[i] = 1.0 / std::sqrt(2.0 * fca_utils::PI_constant()) / sigma * exp(-fca_utils::squareof(1.0 / sigma * (i - std::floor(this->ngauss / 2.0))) / 2);
    }

    /* Initialize to value, even if it is not used */
    this->dvol = -1.;
}

void FcaEntropy2D::setup_bins(const real s[], real* minv, real* maxv, real* dx)
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
    (*maxv) += std::ceil(this->ngauss / 2.0) * (*dx);
    (*minv) -= std::ceil(this->ngauss / 2.0) * (*dx);
    (*dx)    = ((*maxv) - (*minv)) / (1.0 * this->nbins);
}

real FcaEntropy2D::entropy2D(const real s1[], const real s2[])
{
    entropy2D_binning(s1, s2);
    entropy2D_convolution(this->gauss.get(), this->bins.get(), this->p_quant.get());
    return entropy2D_interpolate();
}

} //gmx namespace
