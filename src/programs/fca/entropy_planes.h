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
#ifndef FCA_ENTROPY_PLANES_H
#define FCA_ENTROPY_PLANES_H
#include <math.h>

#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"

#include "utils.hpp"

namespace FCA
{

class Entropy
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

    public:
        Entropy(const int inNdata, const int inNbins)
        {
            this->nbins = inNbins;
            this->ndata = inNdata;

            this->bins.reset(new real[nbins]);
            this->p_quant.reset(new real[nbins]);
            this->phi_quant.reset(new real[nbins]);
            this->wbins.reset(new real[nbins]);
            this->eta.reset(new real[ndata]);
            this->binndx.reset(new int[ndata]);
            this->p.reset(new real[ndata]);
            this->phi.reset(new real[ndata]);

            const real sigma = 5 * (1 + log10(ceil(20000.0 / ndata))); /* fuer 2000 Punkte ist 5 zu niedrig, aber bei 20000 ist es super */
            this->ngauss     = ceil(sigma * 7);                        /* mit *14 ist es genauer, aber so reicht es eigentlich. cumsum(gradS) weicht zwar ein bisschen von S ab, aber hat selben hoch+tief punkte. */

            this->gauss.reset(new real[this->ngauss]);
            this->gaussprime.reset(new real[this->ngauss]);
            this->Fquant.reset(new real[nbins]);
            this->Hgradl.reset(new real[ndata]);

            for (int i = 0; i < this->ngauss; i++)
            {
                this->gauss[i]      = 1.0 / sqrt(2.0 * utils::PI_constant()) / sigma * exp(-utils::squareof(1.0 / sigma * (i - floor(this->ngauss / 2.0))) / 2);
                this->gaussprime[i] = -1.0 / sqrt(2.0 * utils::PI_constant()) / sigma * (1.0 / utils::squareof(sigma) * (i - floor(this->ngauss / 2.0))) * exp(-utils::squareof(1.0 / sigma * (i - floor(this->ngauss / 2.0))) / 2);
            } /* muss dann jeweils mit dx malgenommen werden */
        }

        virtual ~Entropy()
        {
        }

        void entropy_binning(const real s[])
        {
            /*  the data s is binned (number of bins: this->nbins )
               the result is stored in entr
               eta -- interpolation distance to next left bin-edge 0<eta<1  eta=s[i]-binndx[i]
               binndx[i] -- in which bin is s[i]
               bins[j] -- the result of binning.
             */
            real maxv = FCA::utils::limits_min<real>();
            real minv = FCA::utils::limits_max<real>();

            /* find maxima and minima of data */
            for (int i = 0; i < this->ndata; i++)
            {
                maxv = s[i] > maxv ? s[i] : maxv;
                minv = s[i] < minv ? s[i] : minv;
            }

            /* compute bin-width dx */
            real dx = (maxv - minv) / (1.0 * this->nbins);

            /* correct for extra-space at the ends due to gauss-smoothing */
            /* correct to accomodate tails of gaussian-kernel */
            maxv    += ceil(this->ngauss / 2.0) * dx;
            minv    -= ceil(this->ngauss / 2.0) * dx;
            dx       = (maxv - minv) / (1.0 * this->nbins);
            this->dx = dx;

            /* initialize bins with zero */
            for (int i = 0; i < this->nbins; i++)
            {
                this->bins[i]  = 0.0;
                this->wbins[i] = 0.0;
            }

            /* bin the data */
            for (int i = 0; i < this->ndata; i++)
            {
                long int   mhash  = lrint(floor((s[i] - minv) / dx));
                const real eta    = (s[i] - minv) / dx - mhash;
                this->binndx[i]        = mhash;
                this->eta[i]           = 1.0 - eta;
                this->bins[mhash]     += (1.0 - eta);
                this->bins[mhash + 1] += eta;
            }
        }

        void entropy_convolution(const real kernel[], const real bins[], real result[]) const
        {
            /* this computes a convolution of bins with kernel
               could be replaced by fft convolution -- will that be faster ?
               now: 200 bins, ca 35 ngauss 200*35 operations: 7000 (direkt) vs 1035 (fft) but more overhead */
            //const int m      = floor(this->ngauss / 2.0);
            real      invn   = 1.0 / this->ndata;
            for (int i = 0; i < this->nbins; i++)
            {
                result[i]        = 0;
                const int gstart    = std::max(i - this->nbins, 0);
                const int gstopp    = std::min(i + 1, this->ngauss);
                const int shiftbins = std::min(i, this->nbins-1);
                for (int k = gstart; k < gstopp; k++)
                {
                    result[i] += bins[shiftbins - k] * kernel[k] * invn;
                }
                // OK! smoothing verglichen mit MATLAB gaussav.m
            }
        }

        real entropy_interpolate_p()
        {
            real       S            = 0.0;
            const real invdx        = 1.0 / this->dx;
            const int* mhash        = this->binndx.get();
            for (int i = 0; i < this->ndata; i++)
            {
                const real a = this->eta[i];
                const real b = 1 - a;
                this->p[i]                 = (a * this->p_quant[mhash[i]] + b * this->p_quant[mhash[i] + 1]) * invdx;
                S                         += -logf(this->p[i]);
                this->wbins[mhash[i]]     += -a / this->p[i] / this->ndata;
                this->wbins[mhash[i] + 1] += -b / this->p[i] / this->ndata;
            }

            return S * 1.0 / this->ndata;
        }

        void entropy_interpolate_phi()
        {
            const real invdx        = 1.0 / this->dx;
            const real inv2dx       = utils::squareof(invdx);
            const int* mhash        = this->binndx.get();
            for (int i = 0; i < this->ndata; i++)
            {
                const real a = this->eta[i];
                const real b = 1 - a;
                this->phi[i]               = (a * this->phi_quant[mhash[i]] + b * this->phi_quant[mhash[i] + 1]) * inv2dx;
                this->Hgradl[i]            = -this->phi[i] / this->p[i] / this->ndata;
            }

        }

        real entropy_grad(const real gfun[])
        {
            const real invdx        = 1.0 / this->dx;
            const real inv2dx       = utils::squareof(invdx);
            real       gS           = 0;
            for (int i = 0; i < this->ndata; i++)
            {
                const real a = this->eta[i];
                const real b = 1 - a;
                this->Hgradl[i] += (a * this->Fquant[this->binndx[i]] + b * this->Fquant[this->binndx[i] + 1]) * inv2dx;
                gS              += gfun[i] * this->Hgradl[i];
            }
            ;
            return gS;
        }

        /* results:
           S -- resulting 1D entropy
           gS -- nullptr or gradient of entropy
           input :
           s -- data-series of length this->ndata
           g -- nullptr or ...
           entr -- entropy object
         */
        real entropy_1D(real* S, real* gS, const real s[], const real g[])
        {
            entropy_binning(s);
            entropy_convolution(this->gauss.get(), this->bins.get(), this->p_quant.get());
            if (g)
            {
                entropy_convolution(this->gaussprime.get(), this->bins.get(), this->phi_quant.get());
            }
            (*S) = entropy_interpolate_p();
            if (g)
            {
                entropy_interpolate_phi();
                entropy_convolution(this->gaussprime.get(), this->wbins.get(), this->Fquant.get());
                if (gS != nullptr)
                {
                    (*gS) = entropy_grad(g);
                }
            }
            return (*S);
        }
};

class Plane
{
    std::unique_ptr< real[] > proj1; // projections for theta=0;
    std::unique_ptr< real[] > proj2;
    std::unique_ptr< real[] > s1;    // projections for rotation in plane by theta
    std::unique_ptr< real[] > s2;    //s1,s2 for the normal projection
    std::unique_ptr< real[] > g1;    //for gradient
    std::unique_ptr< real[] > g2;
    real theta;
    real S;
    real gS;
    int  ndata;

    public:
        explicit Plane(const int inNdata)
        {
            this->ndata = inNdata;
            this->proj1.reset(new real[ndata]);
            this->proj2.reset(new real[ndata]);
            this->s1.reset(new real[ndata]);
            this->s2.reset(new real[ndata]);
            this->g1.reset(new real[ndata]);
            this->g2.reset(new real[ndata]);
            this->theta = 0;
            this->S     = 0;
            this->gS    = 0;
        }

        virtual ~Plane()
        {
        }

        const real* getS1() const
        {
            return s1.get();
        }

        const real* getS2() const
        {
            return s2.get();
        }

        void init_plane(const std::unique_ptr< real[] > projx[], const real a1[], const real a2[], const int dim)
        {
            for (int i = 0; i < this->ndata; i++)
            {
                this->proj1[i] = 0;
                this->proj2[i] = 0;
                for (int j = 0; j < dim; j++)
                {
                    this->proj1[i] += a1[j] * projx[i][j];
                    this->proj2[i] += a2[j] * projx[i][j];
                }
            }
        }

        void rotate_plane(const real theta)
        {
            /* berechnet projektion auf rotierte achsen */
            /* checked with matlab */
            const real costh = cos(theta);
            const real sinth = sin(theta);
            for (int i = 0; i < this->ndata; i++)
            {
                const real v1c = this->proj1[i] * costh;
                const real v2s = this->proj2[i] * sinth;
                const real v1s = this->proj1[i] * sinth;
                const real v2c = this->proj2[i] * costh;
                this->s1[i]    = v1c + v2s;
                this->s2[i]    = -v1s + v2c;
                this->g1[i]    = -v1s + v2c;
                this->g2[i]    = -v1c - v2s;
            }
        }

        void entropy_plane(real* S, real* gS, Entropy* entr)
        {
            real S1  = 0, S2 = 0;
            real gS1 = 0, gS2 = 0;
            entr->entropy_1D(&S1, &gS1, this->s1.get(), this->g1.get());
            entr->entropy_1D(&S2, &gS2, this->s2.get(), this->g2.get());
            (*S)  = S1 + S2;
            (*gS) = gS1 + gS2;
        }

        real minimize_plane(Entropy* entr)
        {
            const int    SCAN_STEPS = 20;
            const double TOL        = 0.0000001;

            real         S[SCAN_STEPS];
            real         gS[SCAN_STEPS];
            int          k_min = 0;
            real         theta;

            for (int k = 0; k < SCAN_STEPS; k++)
            {
                theta = utils::PI_constant() / SCAN_STEPS / 2 * k;
                rotate_plane(theta);
                entropy_plane(&S[k], &gS[k], entr);
                //  fprintf(stdout,"%10.5f %15.12f %15.12f\n",theta,S[k],gS[k]);
                if (S[k_min] > S[k])
                {
                    k_min = k;
                }
            }

            /* define bracket for brents-method */
            const real bx = utils::PI_constant() / SCAN_STEPS / 2 * k_min;
            const real ax = utils::PI_constant() / SCAN_STEPS / 2 * (k_min - 1);
            const real cx = utils::PI_constant() / SCAN_STEPS / 2 * (k_min + 1);
            //  fprintf(stdout,"initial bracket [%10.5f %10.5f] %10.5f with S=%20.15f\n",ax,cx,bx,S[k_min]);
            const real newS = dbrent_plane(ax, bx, cx, entr, TOL, &theta);
            if (newS > S[k_min])
            {
                fprintf(stdout, "found minimum at [%10.5f] %10.5f Sscan=%10.5f Sold=%10.5f\n", theta, newS, S[k_min], S[0]);
            }

            rotate_plane(theta);
            real newS2;
            entropy_plane(&newS2, &gS[0], entr);
            if ((newS2 - newS) > 1e-6)
            {
                fprintf(stdout, "found minimum at [%10.5f] %20.15f S=%20.15f deltaS %10.5g\n", theta, newS2, newS, newS2 - newS);
            }
            return theta;
            /* bleibt manchmal in lokalen minima stecken. Aber lohnt es sich dann nochmal weiterzusuchen?
               man knnte rund ums gefundene minimum nochmal scannen und dann neue bracket whlen. */
            /*
               S[0]=S[k_min];
               for (k=1,k_min=0;k<SCAN_STEPS;k++) {
               real th=theta+(k-SCAN_STEPS/2)*0.01;
               rotate_plane(plane,th);
               entropy_plane(&S[k],&gS[k],plane,entr);
               fprintf(stdout,"%10.5f %15.12f %15.12f\n",theta,S[k],gS[k]);
               if (S[k_min]>S[k])
               k_min=k;
               };
               if (k_min>0)
               fprintf(stdout,"found new minimum [%10.5f] %10.5f with S=%20.15f\n",theta,theta+(k_min-SCAN_STEPS/2)*0.01,S[k_min]);
             */
        }

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

        real dbrent_plane(const real ax, const real bx, const real cx, Entropy* entr, const real tol, real* xmin)
        {
            const int    ITMAX   = 100;
            const double ZEPS    = 1.0e-10;

            real         d = 1; /* was not initialized , what should it be */
            real         a = (ax < cx ? ax : cx);
            real         b = (ax > cx ? ax : cx);
            real         x = bx;
            real         w = bx;
            real         v = bx;
            rotate_plane(x);

            real fx;
            real dx;
            entropy_plane(&fx, &dx, entr);
            real fw = fx;
            real fv = fx;
            real dw = dx;
            real dv = dx;
            real e  = 0.0;

            for (int iter = 1; iter <= ITMAX; iter++)
            {
                const real xm   = 0.5 * (a + b);
                const real tol1 = tol * fabs(x) + ZEPS;
                const real tol2 = 2.0 * tol1;
                if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
                {
                    *xmin = x;
                    return fx;
                }

                real u;
                if (fabs(e) > tol1)
                {
                    real d1 = 2.0 * (b - a);
                    real d2 = d1;
                    if (dw != dx)
                    {
                        d1 = (w - x) * dx / (dx - dw);
                    }
                    if (dv != dx)
                    {
                        d2 = (v - x) * dx / (dx - dv);
                    }
                    const real u1    = x + d1;
                    const real u2    = x + d2;
                    const int  ok1   = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
                    const int  ok2   = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
                    const real olde  = e;
                    e               = d;
                    if (ok1 || ok2)
                    {
                        if (ok1 && ok2)
                        {
                            d = (fabs(d1) < fabs(d2) ? d1 : d2);
                        }
                        else if (ok1)
                        {
                            d = d1;
                        }
                        else
                        {
                            d = d2;
                        }
                        if (fabs(d) <= fabs(0.5 * olde))
                        {
                            u = x + d;
                            if (u - a < tol2 || b - u < tol2)
                            {
                                d = SIGN(tol1, xm - x);
                            }
                        }
                        else
                        {
                            d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                        }
                    }
                    else
                    {
                        d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                    }
                }
                else
                {
                    d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                }

                real fu;
                real du;
                if (fabs(d) >= tol1)
                {
                    u = x + d;
                    rotate_plane(u);
                    entropy_plane(&fu, &du, entr);
                }
                else
                {
                    u = x + SIGN(tol1, d);
                    rotate_plane(u);
                    entropy_plane(&fu, &du, entr);
                    if (fu > fx)
                    {
                        *xmin = x;
                        return fx;
                    }
                }

                if (fu <= fx)
                {
                    if (u >= x)
                    {
                        a = x;
                    }
                    else
                    {
                        b = x;
                    }
                    MOV3(v, fv, dv, w, fw, dw);
                    MOV3(w, fw, dw, x, fx, dx);
                    MOV3(x, fx, dx, u, fu, du);
                }
                else
                {
                    if (u < x)
                    {
                        a = u;
                    }
                    else
                    {
                        b = u;
                    }
                    if (fu <= fw || w == x)
                    {
                        MOV3(v, fv, dv, w, fw, dw);
                        MOV3(w, fw, dw, u, fu, du);
                    }
                    else if (fu < fv || v == x || v == w)
                    {
                        MOV3(v, fv, dv, u, fu, du);
                    }
                }
            }
            fprintf(stdout, "FATAL: Too many iterations in routine DBRENT");
            exit(1);
            return FCA::utils::limits_max<real>();
        }
};
}

#endif
