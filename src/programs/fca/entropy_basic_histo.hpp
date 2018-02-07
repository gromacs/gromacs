#ifndef ENTROPY_BASIC_HISTO_HPP
#define ENTROPY_BASIC_HISTO_HPP

/*
 * $Id: g_project.c,v 1.5 2004/11/12 17:17:36 olange Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to  consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */
//static char *SRCID_g_anaeig_c = "$Id: g_project.c,v 1.5 2004/11/12 17:17:36 olange Exp $";
#include <cmath>
#include <cassert>
#include <memory>
#include "gromacs/utility/futil.h"
#include "fca_minimizing.h"

namespace FCA {

namespace EntropyHisto {

/* results:
   S -- resulting 1D entropy
   input :
   s -- data-series of length entr->ndata
   entr -- entropy object
*/
struct BasicEntropy {
    real log_kappa1; //proportionality for binwidth
    real log_kappa2; //prop. for 2D bandwith
    int ndata;
};


inline real entropy_1D_basic_histo(real* pS, const real s[], const BasicEntropy& entr) {
    assert(pS);
    real maxv         = FCA::utils::limits_min<real>();
    real minv         = FCA::utils::limits_max<real>();
    real avg          = 0;
    real var          = 0;
    const double invn = 1.0 / entr.ndata;
    /* find maxima and minima of data */
    for(int i = 0; i < entr.ndata; i++) {
        maxv = FCA::utils::maxv(s[i], maxv);
        //fprintf(stdout,"maxv %f minv %f s %f\n",maxv,minv,s[i]);
        minv = FCA::utils::minv(s[i], minv);
        avg += s[i] * invn;
        var += s[i] * s[i] * invn;
    };

    const real sig      = var - avg * avg;
    const real binwidth = exp(entr.log_kappa1) * sqrt(sig);
    const int nbins     = lrint(ceil((maxv - minv) / binwidth)) + 1;
    real norm           = binwidth;

    std::unique_ptr< real[] > bins(new real[nbins]());
    //bin data
    for(int i = 0; i < entr.ndata; i++) {
        const int ibin = lrint(floor((s[i] - minv) / binwidth));
        bins[ibin] += invn;
    }

    //integrate bins
    real S = 0;
    for(int i = 0; i < nbins; i++) {
        if(bins[i] > 0.0) {
            S -= bins[i] * log(bins[i] / norm);
        }
    }

    (*pS) = S;
    return (*pS);
}

inline real entropy_2D_basic_histo(const real sx[], const real sy[], const BasicEntropy& entr) {
    real x_max = FCA::utils::limits_min<real>();
    real x_min = FCA::utils::limits_max<real>();
    real y_max = FCA::utils::limits_min<real>();
    real y_min = FCA::utils::limits_max<real>();

    real xavg = 0;
    real xvar = 0;
    real yavg = 0;
    real yvar = 0;

    const real invn = 1.0 / entr.ndata;
    /* find maxima and minima of data */
    for(int i = 0; i < entr.ndata; i++) {
        x_max = FCA::utils::maxv(sx[i], x_max);
        //fprintf(sxtderr,"x_max %f x_min %f sx %f\n",x_max,x_min,sx[i]);
        x_min = FCA::utils::maxv(sx[i], x_min);

        y_max = FCA::utils::maxv(sy[i], y_max);
        //fprintf(sytderr,"y_max %f y_min %f sy %f\n",y_max,y_min,sy[i]);
        y_min = FCA::utils::minv(sy[i], y_min);

        xavg += sx[i] * invn;
        xvar += sx[i] * sx[i] * invn;

        yavg += sy[i] * invn;
        yvar += sy[i] * sy[i] * invn;
    };

    const real x_sig      = sqrt(xvar - xavg * xavg);
    const real x_binwidth = exp(entr.log_kappa2) * x_sig;
    const int x_nbins     = lrint(ceil((x_max - x_min) / x_binwidth)) + 1;

    const real y_sig      = sqrt(yvar - yavg * yavg);
    const real y_binwidth = exp(entr.log_kappa2) * y_sig;
    const int y_nbins     = lrint(ceil((y_max - y_min) / y_binwidth)) + 1;

    //  fprintf( stdout, "bins : %d %d width %5.3f %5.3f sig %5.3f %5.3f ", x_nbins, y_nbins, x_binwidth, y_binwidth, x_sig, y_sig );
    //  real *bins=nullptr;
    //  snew( bins, x_nbins * y_nbins );
    std::unique_ptr< real[] > bins(new real[x_nbins * y_nbins]());

    //bin data
    for(int i = 0; i < entr.ndata; i++) {
        const int ix = lrint(floor((sx[i] - x_min) / x_binwidth));
        const int iy = lrint(floor((sy[i] - y_min) / y_binwidth));
        bins[ix + iy * x_nbins] += invn;
    }
    //integrate bins
    real S          = 0;
    const real norm = x_binwidth * y_binwidth;
    for(int i = 0; i < x_nbins; i++) {
        for(int j = 0; j < y_nbins; j++) {
            const real b = bins[i + j * x_nbins];
            if(b > 0) {
                S += -b * log(b / norm);
            }
        }
    }
    return S;
}

inline void Basic_compute_MI_matrix(FcaMaster* fca, const real log_kappa1, const real log_kappa2) {
    BasicEntropy entr;
    entr.ndata      = fca->getNFrames();
    entr.log_kappa1 = log_kappa1;
    entr.log_kappa2 = log_kappa2;

    std::unique_ptr< real[] >* icx = fca->getIcx();
    const int dim                  = fca->getDim();
    real* S1                       = fca->getS1();

    {
        FILE* S1D_file = gmx_ffopen("entropy_basic_1D.dat", "w");
        for(int i = 0; i < dim; i++) {
            entropy_1D_basic_histo(&S1[i], icx[i].get(), entr);
            fprintf(S1D_file, "%5d %10.7f\n", i, S1[i]);
        }
        fclose(S1D_file);
    }

    double* MIs = fca->getMI();
    std::unique_ptr<double[]> S2D(new double[dim * dim]());
    for(int i = 0; i < dim; i++) {
        for(int j = i; j < dim; j++) {
            if(i != j) {
                real S2          = entropy_2D_basic_histo(icx[i].get(), icx[j].get(), entr);
                real S1s         = S1[i] + S1[j];
                MIs[i + dim * j] = S1s - S2;
                MIs[j + dim * i] = S1s - S2;
                fprintf(stdout, "SINGLE_MODE %2d %2d %5.3f %5.3f %5.3f\n", i, j, S2, S1[i], S1[j]);
                S2D[j + dim * i] = S2;
                S2D[i + dim * j] = S2;
            } else {
                MIs[i + dim * j] = 0.0;
                S2D[i + dim * i] = 0.0;
            }
        }
    }

    {
        FILE* S2D_file = gmx_ffopen("entropy_basic_2D.dat", "w");
        utils::dump_sqrmatrix(S2D_file, S2D.get(), dim);
        fclose(S2D_file);
    }
}

}

}


#endif
