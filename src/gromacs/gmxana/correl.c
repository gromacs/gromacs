/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gromacs/fft/fft.h"
#include "gromacs/utility/smalloc.h"
#include "correl.h"

#define SWAP(a, b) tempr = (a); (a) = (b); (b) = tempr

void four1(real data[], int nn, int isign)
{
    int    n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    real   tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m)
        {
            j  -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax)
    {
        istep = 2*mmax;
        theta = 6.28318530717959/(isign*mmax);
        wtemp = sin(0.5*theta);
        wpr   = -2.0*wtemp*wtemp;
        wpi   = sin(theta);
        wr    = 1.0;
        wi    = 0.0;
        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j          = i+mmax;
                tempr      = wr*data[j]-wi*data[j+1];
                tempi      = wr*data[j+1]+wi*data[j];
                data[j]    = data[i]-tempr;
                data[j+1]  = data[i+1]-tempi;
                data[i]   += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr)*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
        }
        mmax = istep;
    }
}

#undef SWAP

static void realft(real data[], int n, int isign)
{
    int    i, i1, i2, i3, i4, n2p3;
    real   c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    double wr, wi, wpr, wpi, wtemp, theta;

    theta = 3.141592653589793/(double) n;
    if (isign == 1)
    {
        c2 = -0.5;
        four1(data, n, 1);
    }
    else
    {
        c2    = 0.5;
        theta = -theta;
    }
    wtemp = sin(0.5*theta);
    wpr   = -2.0*wtemp*wtemp;
    wpi   = sin(theta);
    wr    = 1.0+wpr;
    wi    = wpi;
    n2p3  = 2*n+3;
    for (i = 2; i <= n/2; i++)
    {
        i4       = 1+(i3 = n2p3-(i2 = 1+(i1 = i+i-1)));
        h1r      = c1*(data[i1]+data[i3]);
        h1i      = c1*(data[i2]-data[i4]);
        h2r      = -c2*(data[i2]+data[i4]);
        h2i      = c2*(data[i1]-data[i3]);
        data[i1] = h1r+wr*h2r-wi*h2i;
        data[i2] = h1i+wr*h2i+wi*h2r;
        data[i3] = h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr       = (wtemp = wr)*wpr-wi*wpi+wr;
        wi       = wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1)
    {
        data[1] = (h1r = data[1])+data[2];
        data[2] = h1r-data[2];
    }
    else
    {
        data[1] = c1*((h1r = data[1])+data[2]);
        data[2] = c1*(h1r-data[2]);
        four1(data, n, -1);
    }
}

static void twofft(real data1[], real data2[], real fft1[], real fft2[], int n)
{
    int  nn3, nn2, jj, j;
    real rep, rem, aip, aim;

    nn3 = 1+(nn2 = 2+n+n);
    for (j = 1, jj = 2; j <= n; j++, jj += 2)
    {
        fft1[jj-1] = data1[j];
        fft1[jj]   = data2[j];
    }
    four1(fft1, n, 1);
    fft2[1] = fft1[2];
    fft1[2] = fft2[2] = 0.0;
    for (j = 3; j <= n+1; j += 2)
    {
        rep         = 0.5*(fft1[j]+fft1[nn2-j]);
        rem         = 0.5*(fft1[j]-fft1[nn2-j]);
        aip         = 0.5*(fft1[j+1]+fft1[nn3-j]);
        aim         = 0.5*(fft1[j+1]-fft1[nn3-j]);
        fft1[j]     = rep;
        fft1[j+1]   = aim;
        fft1[nn2-j] = rep;
        fft1[nn3-j] = -aim;
        fft2[j]     = aip;
        fft2[j+1]   = -rem;
        fft2[nn2-j] = aip;
        fft2[nn3-j] = rem;
    }
}

void correl(real data1[], real data2[], int n, real ans[])
{
    int     no2, i;
    real    dum, *fft;

    snew(fft, 2*n+1);
    twofft(data1, data2, fft, ans, n);
    no2 = n/2;
    for (i = 2; i <= n+2; i += 2)
    {
        dum      = ans[i-1];
        ans[i-1] = (fft[i-1]*dum+fft[i]*ans[i])/no2;
        ans[i]   = (fft[i]*dum-fft[i-1]*ans[i])/no2;
    }
    ans[2] = ans[n+1];
    realft(ans, no2, -1);
    sfree(fft);
}
