/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_FFT_FFTPACK

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


#include "gmx_fft.h"
#include "gmx_fatal.h"


/** Contents of the FFTPACK fft datatype. 
 *
 *  FFTPACK only does 1d transforms, so we use a pointers to another fft for 
 *  the transform in the next dimension.
 * Thus, a 3d-structure contains a pointer to a 2d one, which in turns contains
 * a pointer to a 1d. The 1d structure has next==NULL.
 */
struct gmx_fft
{
    int            ndim;     /**< Dimensions, including our subdimensions.  */
    int            n;        /**< Number of points in this dimension.       */
    int            ifac[15]; /**< 15 bytes needed for cfft and rfft         */
    struct gmx_fft *next;    /**< Pointer to next dimension, or NULL.       */
    real *         work;     /**< 1st 4n reserved for cfft, 1st 2n for rfft */
};

#include <math.h>
#include <stdio.h>



static void 
fftpack_passf2(int         ido, 
               int         l1, 
               real  cc[],
               real  ch[],
               real  wa1[],
               int         isign)
{
    int i, k, ah, ac;
    real ti2, tr2;
    
    if (ido <= 2)
    {
        for (k=0; k<l1; k++) 
        {
            ah = k*ido;
            ac = 2*k*ido;
            ch[ah]              = cc[ac]   + cc[ac + ido];
            ch[ah + ido*l1]     = cc[ac]   - cc[ac + ido];
            ch[ah+1]            = cc[ac+1] + cc[ac + ido + 1];
            ch[ah + ido*l1 + 1] = cc[ac+1] - cc[ac + ido + 1];
        }
    } 
    else
    {
        for (k=0; k<l1; k++) 
        {
            for (i=0; i<ido-1; i+=2) 
            {
                ah              = i + k*ido;
                ac              = i + 2*k*ido;
                ch[ah]          = cc[ac] + cc[ac + ido];
                tr2             = cc[ac] - cc[ac + ido];
                ch[ah+1]        = cc[ac+1] + cc[ac + 1 + ido];
                ti2             = cc[ac+1] - cc[ac + 1 + ido];
                ch[ah+l1*ido+1] = wa1[i]*ti2 + isign*wa1[i+1]*tr2;
                ch[ah+l1*ido]   = wa1[i]*tr2 - isign*wa1[i+1]*ti2;
            }
        }
    }
} 



static void 
fftpack_passf3(int         ido,
               int         l1, 
               real  cc[],
               real  ch[],
               real  wa1[], 
               real  wa2[],
               int         isign)
{
    const real taur = -0.5;
    const real taui = 0.866025403784439;
    
    int i, k, ac, ah;
    real ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    
    if (ido == 2) 
    {
        for (k=1; k<=l1; k++) 
        {
            ac = (3*k - 2)*ido;
            tr2 = cc[ac] + cc[ac + ido];
            cr2 = cc[ac - ido] + taur*tr2;
            ah = (k - 1)*ido;
            ch[ah] = cc[ac - ido] + tr2;
            
            ti2 = cc[ac + 1] + cc[ac + ido + 1];
            ci2 = cc[ac - ido + 1] + taur*ti2;
            ch[ah + 1] = cc[ac - ido + 1] + ti2;
            
            cr3 = isign*taui*(cc[ac] - cc[ac + ido]);
            ci3 = isign*taui*(cc[ac + 1] - cc[ac + ido + 1]);
            ch[ah + l1*ido] = cr2 - ci3;
            ch[ah + 2*l1*ido] = cr2 + ci3;
            ch[ah + l1*ido + 1] = ci2 + cr3;
            ch[ah + 2*l1*ido + 1] = ci2 - cr3;
        }
    } 
    else
    {
        for (k=1; k<=l1; k++) 
        {
            for (i=0; i<ido-1; i+=2)
            {
                ac = i + (3*k - 2)*ido;
                tr2 = cc[ac] + cc[ac + ido];
                cr2 = cc[ac - ido] + taur*tr2;
                ah = i + (k-1)*ido;
                ch[ah] = cc[ac - ido] + tr2;
                ti2 = cc[ac + 1] + cc[ac + ido + 1];
                ci2 = cc[ac - ido + 1] + taur*ti2;
                ch[ah + 1] = cc[ac - ido + 1] + ti2;
                cr3 = isign*taui*(cc[ac] - cc[ac + ido]);
                ci3 = isign*taui*(cc[ac + 1] - cc[ac + ido + 1]);
                dr2 = cr2 - ci3;
                dr3 = cr2 + ci3;
                di2 = ci2 + cr3;
                di3 = ci2 - cr3;
                ch[ah + l1*ido + 1] = wa1[i]*di2 + isign*wa1[i+1]*dr2;
                ch[ah + l1*ido] = wa1[i]*dr2 - isign*wa1[i+1]*di2;
                ch[ah + 2*l1*ido + 1] = wa2[i]*di3 + isign*wa2[i+1]*dr3;
                ch[ah + 2*l1*ido] = wa2[i]*dr3 - isign*wa2[i+1]*di3;
            }
        }
    }
} 


static void 
fftpack_passf4(int          ido, 
               int          l1,
               real   cc[], 
               real   ch[],
               real   wa1[],
               real   wa2[], 
               real   wa3[],
               int          isign)
{
    int i, k, ac, ah;
    real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    
    if (ido == 2) 
    {
        for (k=0; k<l1; k++)
        {
            ac = 4*k*ido + 1;
            ti1 = cc[ac] - cc[ac + 2*ido];
            ti2 = cc[ac] + cc[ac + 2*ido];
            tr4 = cc[ac + 3*ido] - cc[ac + ido];
            ti3 = cc[ac + ido] + cc[ac + 3*ido];
            tr1 = cc[ac - 1] - cc[ac + 2*ido - 1];
            tr2 = cc[ac - 1] + cc[ac + 2*ido - 1];
            ti4 = cc[ac + ido - 1] - cc[ac + 3*ido - 1];
            tr3 = cc[ac + ido - 1] + cc[ac + 3*ido - 1];
            ah = k*ido;
            ch[ah] = tr2 + tr3;
            ch[ah + 2*l1*ido] = tr2 - tr3;
            ch[ah + 1] = ti2 + ti3;
            ch[ah + 2*l1*ido + 1] = ti2 - ti3;
            ch[ah + l1*ido] = tr1 + isign*tr4;
            ch[ah + 3*l1*ido] = tr1 - isign*tr4;
            ch[ah + l1*ido + 1] = ti1 + isign*ti4;
            ch[ah + 3*l1*ido + 1] = ti1 - isign*ti4;
        }
    } 
    else
    {
        for (k=0; k<l1; k++)
        {
            for (i=0; i<ido-1; i+=2)
            {
                ac = i + 1 + 4*k*ido;
                ti1 = cc[ac] - cc[ac + 2*ido];
                ti2 = cc[ac] + cc[ac + 2*ido];
                ti3 = cc[ac + ido] + cc[ac + 3*ido];
                tr4 = cc[ac + 3*ido] - cc[ac + ido];
                tr1 = cc[ac - 1] - cc[ac + 2*ido - 1];
                tr2 = cc[ac - 1] + cc[ac + 2*ido - 1];
                ti4 = cc[ac + ido - 1] - cc[ac + 3*ido - 1];
                tr3 = cc[ac + ido - 1] + cc[ac + 3*ido - 1];
                ah = i + k*ido;
                ch[ah] = tr2 + tr3;
                cr3 = tr2 - tr3;
                ch[ah + 1] = ti2 + ti3;
                ci3 = ti2 - ti3;
                cr2 = tr1 + isign*tr4;
                cr4 = tr1 - isign*tr4;
                ci2 = ti1 + isign*ti4;
                ci4 = ti1 - isign*ti4;
                ch[ah + l1*ido] = wa1[i]*cr2 - isign*wa1[i + 1]*ci2;
                ch[ah + l1*ido + 1] = wa1[i]*ci2 + isign*wa1[i + 1]*cr2;
                ch[ah + 2*l1*ido] = wa2[i]*cr3 - isign*wa2[i + 1]*ci3;
                ch[ah + 2*l1*ido + 1] = wa2[i]*ci3 + isign*wa2[i + 1]*cr3;
                ch[ah + 3*l1*ido] = wa3[i]*cr4 -isign*wa3[i + 1]*ci4;
                ch[ah + 3*l1*ido + 1] = wa3[i]*ci4 + isign*wa3[i + 1]*cr4;
            }
        }
    }
} 


static void 
fftpack_passf5(int          ido,
               int          l1, 
               real   cc[],
               real   ch[],
               real   wa1[], 
               real   wa2[],
               real   wa3[], 
               real   wa4[],
               int          isign)
{
    const real tr11 = 0.309016994374947;
    const real ti11 = 0.951056516295154;
    const real tr12 = -0.809016994374947;
    const real ti12 = 0.587785252292473;
    
    int i, k, ac, ah;
    real ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
        ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    
    if (ido == 2) 
    {
        for (k = 1; k <= l1; ++k) 
        {
            ac = (5*k - 4)*ido + 1;
            ti5 = cc[ac] - cc[ac + 3*ido];
            ti2 = cc[ac] + cc[ac + 3*ido];
            ti4 = cc[ac + ido] - cc[ac + 2*ido];
            ti3 = cc[ac + ido] + cc[ac + 2*ido];
            tr5 = cc[ac - 1] - cc[ac + 3*ido - 1];
            tr2 = cc[ac - 1] + cc[ac + 3*ido - 1];
            tr4 = cc[ac + ido - 1] - cc[ac + 2*ido - 1];
            tr3 = cc[ac + ido - 1] + cc[ac + 2*ido - 1];
            ah = (k - 1)*ido;
            ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
            ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
            cr2 = cc[ac - ido - 1] + tr11*tr2 + tr12*tr3;
            ci2 = cc[ac - ido] + tr11*ti2 + tr12*ti3;
            cr3 = cc[ac - ido - 1] + tr12*tr2 + tr11*tr3;
            ci3 = cc[ac - ido] + tr12*ti2 + tr11*ti3;
            cr5 = isign*(ti11*tr5 + ti12*tr4);
            ci5 = isign*(ti11*ti5 + ti12*ti4);
            cr4 = isign*(ti12*tr5 - ti11*tr4);
            ci4 = isign*(ti12*ti5 - ti11*ti4);
            ch[ah + l1*ido] = cr2 - ci5;
            ch[ah + 4*l1*ido] = cr2 + ci5;
            ch[ah + l1*ido + 1] = ci2 + cr5;
            ch[ah + 2*l1*ido + 1] = ci3 + cr4;
            ch[ah + 2*l1*ido] = cr3 - ci4;
            ch[ah + 3*l1*ido] = cr3 + ci4;
            ch[ah + 3*l1*ido + 1] = ci3 - cr4;
            ch[ah + 4*l1*ido + 1] = ci2 - cr5;
        }
    } 
    else
    {
        for (k=1; k<=l1; k++) 
        {
            for (i=0; i<ido-1; i+=2) 
            {
                ac = i + 1 + (k*5 - 4)*ido;
                ti5 = cc[ac] - cc[ac + 3*ido];
                ti2 = cc[ac] + cc[ac + 3*ido];
                ti4 = cc[ac + ido] - cc[ac + 2*ido];
                ti3 = cc[ac + ido] + cc[ac + 2*ido];
                tr5 = cc[ac - 1] - cc[ac + 3*ido - 1];
                tr2 = cc[ac - 1] + cc[ac + 3*ido - 1];
                tr4 = cc[ac + ido - 1] - cc[ac + 2*ido - 1];
                tr3 = cc[ac + ido - 1] + cc[ac + 2*ido - 1];
                ah = i + (k - 1)*ido;
                ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
                ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
                cr2 = cc[ac - ido - 1] + tr11*tr2 + tr12*tr3;
                ci2 = cc[ac - ido] + tr11*ti2 + tr12*ti3;
                cr3 = cc[ac - ido - 1] + tr12*tr2 + tr11*tr3;
                ci3 = cc[ac - ido] + tr12*ti2 + tr11*ti3;
                cr5 = isign*(ti11*tr5 + ti12*tr4);
                ci5 = isign*(ti11*ti5 + ti12*ti4);
                cr4 = isign*(ti12*tr5 - ti11*tr4);
                ci4 = isign*(ti12*ti5 - ti11*ti4);
                dr3 = cr3 - ci4;
                dr4 = cr3 + ci4;
                di3 = ci3 + cr4;
                di4 = ci3 - cr4;
                dr5 = cr2 + ci5;
                dr2 = cr2 - ci5;
                di5 = ci2 - cr5;
                di2 = ci2 + cr5;
                ch[ah + l1*ido] = wa1[i]*dr2 - isign*wa1[i+1]*di2;
                ch[ah + l1*ido + 1] = wa1[i]*di2 + isign*wa1[i+1]*dr2;
                ch[ah + 2*l1*ido] = wa2[i]*dr3 - isign*wa2[i+1]*di3;
                ch[ah + 2*l1*ido + 1] = wa2[i]*di3 + isign*wa2[i+1]*dr3;
                ch[ah + 3*l1*ido] = wa3[i]*dr4 - isign*wa3[i+1]*di4;
                ch[ah + 3*l1*ido + 1] = wa3[i]*di4 + isign*wa3[i+1]*dr4;
                ch[ah + 4*l1*ido] = wa4[i]*dr5 - isign*wa4[i+1]*di5;
                ch[ah + 4*l1*ido + 1] = wa4[i]*di5 + isign*wa4[i+1]*dr5;
            }
        }
    }
} 


static void 
fftpack_passf(int *        nac, 
              int          ido,
              int          ip, 
              int          l1,
              int          idl1,
              real   cc[],
              real   ch[],
              real   wa[],
              int          isign)
{
    int idij, idlj, idot, ipph, i, j, k, l, jc, lc, ik, nt, idj, idl, inc,idp;
    real wai, war;
    
    idot = ido / 2;
    nt = ip*idl1;
    ipph = (ip + 1) / 2;
    idp = ip*ido;
    if (ido >= l1) 
    {
        for (j=1; j<ipph; j++)
        {
            jc = ip - j;
            for (k=0; k<l1; k++) 
            {
                for (i=0; i<ido; i++) 
                {
                    ch[i + (k + j*l1)*ido]  = cc[i + (j + k*ip)*ido] + cc[i + (jc + k*ip)*ido];
                    ch[i + (k + jc*l1)*ido] = cc[i + (j + k*ip)*ido] - cc[i + (jc + k*ip)*ido];
                }
            }
        }
        for (k=0; k<l1; k++)
            for (i=0; i<ido; i++)
                ch[i + k*ido] = cc[i + k*ip*ido];
    } 
    else
    {
        for (j=1; j<ipph; j++) 
        {
            jc = ip - j;
            for (i=0; i<ido; i++) 
            {
                for (k=0; k<l1; k++) 
                {
                    ch[i + (k + j*l1)*ido] =  cc[i + (j + k*ip)*ido] + cc[i + (jc + k*ip)*ido];
                    ch[i + (k + jc*l1)*ido] = cc[i + (j + k*ip)*ido] - cc[i + (jc + k*ip)*ido];
                }
            }
        }
        for (i=0; i<ido; i++)
            for (k=0; k<l1; k++)
                ch[i + k*ido] = cc[i + k*ip*ido];
    }
    
    idl = 2 - ido;
    inc = 0;
    for (l=1; l<ipph; l++) 
    {
        lc = ip - l;
        idl += ido;
        for (ik=0; ik<idl1; ik++)
        {
            cc[ik + l*idl1] = ch[ik] + wa[idl - 2]*ch[ik + idl1];
            cc[ik + lc*idl1] = isign*wa[idl-1]*ch[ik + (ip-1)*idl1];
        }
        idlj = idl;
        inc += ido;
        for (j=2; j<ipph; j++)
        {
            jc = ip - j;
            idlj += inc;
            if (idlj > idp) idlj -= idp;
            war = wa[idlj - 2];
            wai = wa[idlj-1];
            for (ik=0; ik<idl1; ik++)
            {
                cc[ik + l*idl1] += war*ch[ik + j*idl1];
                cc[ik + lc*idl1] += isign*wai*ch[ik + jc*idl1];
            }
        }
    }
    for (j=1; j<ipph; j++)
        for (ik=0; ik<idl1; ik++)
            ch[ik] += ch[ik + j*idl1];
    for (j=1; j<ipph; j++) 
    {
        jc = ip - j;
        for (ik=1; ik<idl1; ik+=2) 
        {
            ch[ik - 1 + j*idl1] = cc[ik - 1 + j*idl1] - cc[ik + jc*idl1];
            ch[ik - 1 + jc*idl1] = cc[ik - 1 + j*idl1] + cc[ik + jc*idl1];
            ch[ik + j*idl1] = cc[ik + j*idl1] + cc[ik - 1 + jc*idl1];
            ch[ik + jc*idl1] = cc[ik + j*idl1] - cc[ik - 1 + jc*idl1];
        }
    }
    *nac = 1;
    if (ido == 2) 
        return;
    *nac = 0;
    for (ik=0; ik<idl1; ik++)
    {
        cc[ik] = ch[ik];
    }
    for (j=1; j<ip; j++)
    {
        for (k=0; k<l1; k++) 
        {
            cc[(k + j*l1)*ido + 0] = ch[(k + j*l1)*ido + 0];
            cc[(k + j*l1)*ido + 1] = ch[(k + j*l1)*ido + 1];
        }
    }
    if (idot <= l1) 
    {
        idij = 0;
        for (j=1; j<ip; j++)
        {
            idij += 2;
            for (i=3; i<ido; i+=2) 
            {
                idij += 2;
                for (k=0; k<l1; k++)
                {
                    cc[i - 1 + (k + j*l1)*ido] =
                    wa[idij - 2]*ch[i - 1 + (k + j*l1)*ido] -
                    isign*wa[idij-1]*ch[i + (k + j*l1)*ido];
                    cc[i + (k + j*l1)*ido] =
                        wa[idij - 2]*ch[i + (k + j*l1)*ido] +
                        isign*wa[idij-1]*ch[i - 1 + (k + j*l1)*ido];
                }
            }
        }
    }
    else
    {
        idj = 2 - ido;
        for (j=1; j<ip; j++) 
        {
            idj += ido;
            for (k = 0; k < l1; k++) 
            {
                idij = idj;
                for (i=3; i<ido; i+=2)
                {
                    idij += 2;
                    cc[i - 1 + (k + j*l1)*ido] =
                        wa[idij - 2]*ch[i - 1 + (k + j*l1)*ido] -
                        isign*wa[idij-1]*ch[i + (k + j*l1)*ido];
                    cc[i + (k + j*l1)*ido] =
                        wa[idij - 2]*ch[i + (k + j*l1)*ido] +
                        isign*wa[idij-1]*ch[i - 1 + (k + j*l1)*ido];
                }
            }
        }
    }
} 



static void 
fftpack_radf2(int          ido,
              int          l1, 
              real   cc[], 
              real   ch[], 
              real   wa1[])
{
    int i, k, ic;
    real ti2, tr2;
    for (k=0; k<l1; k++) 
    {
        ch[2*k*ido] = cc[k*ido] + cc[(k + l1)*ido];
        ch[(2*k+1)*ido + ido-1] = cc[k*ido] - cc[(k + l1)*ido];
    }
    if (ido < 2) 
        return;
    if (ido != 2) 
    {
        for (k=0; k<l1; k++) 
        {
            for (i=2; i<ido; i+=2) 
            {
                ic = ido - i;
                tr2 = wa1[i - 2]*cc[i-1 + (k + l1)*ido] + wa1[i - 1]*cc[i + (k + l1)*ido];
                ti2 = wa1[i - 2]*cc[i + (k + l1)*ido] - wa1[i - 1]*cc[i-1 + (k + l1)*ido];
                ch[i + 2*k*ido] = cc[i + k*ido] + ti2;
                ch[ic + (2*k+1)*ido] = ti2 - cc[i + k*ido];
                ch[i - 1 + 2*k*ido] = cc[i - 1 + k*ido] + tr2;
                ch[ic - 1 + (2*k+1)*ido] = cc[i - 1 + k*ido] - tr2;
            }
        }
        if (ido % 2 == 1) 
            return;
    }
    for (k=0; k<l1; k++)
    {
        ch[(2*k+1)*ido] = -cc[ido-1 + (k + l1)*ido];
        ch[ido-1 + 2*k*ido] = cc[ido-1 + k*ido];
    }
}


static void 
fftpack_radb2(int          ido, 
              int          l1, 
              real   cc[],
              real   ch[], 
              real   wa1[])
{
    int i, k, ic;
    real ti2, tr2;
    for (k=0; k<l1; k++) 
    {
        ch[k*ido] = cc[2*k*ido] + cc[ido-1 + (2*k+1)*ido];
        ch[(k + l1)*ido] = cc[2*k*ido] - cc[ido-1 + (2*k+1)*ido];
    }
    if (ido < 2) 
        return;
    if (ido != 2) 
    {
        for (k = 0; k < l1; ++k)
        {
            for (i = 2; i < ido; i += 2) 
            {
                ic = ido - i;
                ch[i-1 + k*ido] = cc[i-1 + 2*k*ido] + cc[ic-1 + (2*k+1)*ido];
                tr2 = cc[i-1 + 2*k*ido] - cc[ic-1 + (2*k+1)*ido];
                ch[i + k*ido] = cc[i + 2*k*ido] - cc[ic + (2*k+1)*ido];
                ti2 = cc[i + (2*k)*ido] + cc[ic + (2*k+1)*ido];
                ch[i-1 + (k + l1)*ido] = wa1[i - 2]*tr2 - wa1[i - 1]*ti2;
                ch[i + (k + l1)*ido] = wa1[i - 2]*ti2 + wa1[i - 1]*tr2;
            }
        }
        if (ido % 2 == 1) 
            return;
    }
    for (k = 0; k < l1; k++) 
    {
        ch[ido-1 + k*ido] = 2*cc[ido-1 + 2*k*ido];
        ch[ido-1 + (k + l1)*ido] = -2*cc[(2*k+1)*ido];
    }
} 


static void 
fftpack_radf3(int          ido, 
              int          l1,
              real   cc[], 
              real   ch[],
              real   wa1[], 
              real   wa2[])
{
    const real taur = -0.5;
    const real taui = 0.866025403784439;
    int i, k, ic;
    real ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    
    for (k=0; k<l1; k++)
    {
        cr2 = cc[(k + l1)*ido] + cc[(k + 2*l1)*ido];
        ch[3*k*ido] = cc[k*ido] + cr2;
        ch[(3*k+2)*ido] = taui*(cc[(k + l1*2)*ido] - cc[(k + l1)*ido]);
        ch[ido-1 + (3*k + 1)*ido] = cc[k*ido] + taur*cr2;
    }
    if (ido == 1) 
        return;
    for (k=0; k<l1; k++) 
    {
        for (i=2; i<ido; i+=2) 
        {
            ic = ido - i;
            dr2 = wa1[i - 2]*cc[i - 1 + (k + l1)*ido] +wa1[i - 1]*cc[i + (k + l1)*ido];
            di2 = wa1[i - 2]*cc[i + (k + l1)*ido] - wa1[i - 1]*cc[i - 1 + (k + l1)*ido];
            dr3 = wa2[i - 2]*cc[i - 1 + (k + l1*2)*ido] + wa2[i - 1]*cc[i + (k + l1*2)*ido];
            di3 = wa2[i - 2]*cc[i + (k + l1*2)*ido] - wa2[i - 1]*cc[i - 1 + (k + l1*2)*ido];
            cr2 = dr2 + dr3;
            ci2 = di2 + di3;
            ch[i - 1 + 3*k*ido] = cc[i - 1 + k*ido] + cr2;
            ch[i + 3*k*ido] = cc[i + k*ido] + ci2;
            tr2 = cc[i - 1 + k*ido] + taur*cr2;
            ti2 = cc[i + k*ido] + taur*ci2;
            tr3 = taui*(di2 - di3);
            ti3 = taui*(dr3 - dr2);
            ch[i - 1 + (3*k + 2)*ido] = tr2 + tr3;
            ch[ic - 1 + (3*k + 1)*ido] = tr2 - tr3;
            ch[i + (3*k + 2)*ido] = ti2 + ti3;
            ch[ic + (3*k + 1)*ido] = ti3 - ti2;
        }
    }
} 


static void 
fftpack_radb3(int          ido, 
              int          l1, 
              real   cc[], 
              real   ch[],
              real   wa1[],
              real   wa2[])
{
    const real taur = -0.5;
    const real taui = 0.866025403784439;
    int i, k, ic;
    real ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    
    for (k=0; k<l1; k++) 
    {
        tr2 = 2*cc[ido-1 + (3*k + 1)*ido];
        cr2 = cc[3*k*ido] + taur*tr2;
        ch[k*ido] = cc[3*k*ido] + tr2;
        ci3 = 2*taui*cc[(3*k + 2)*ido];
        ch[(k + l1)*ido] = cr2 - ci3;
        ch[(k + 2*l1)*ido] = cr2 + ci3;
    }
    if (ido == 1) 
        return;
    
    for (k=0; k<l1; k++)
    {
        for (i=2; i<ido; i+=2) 
        {
            ic = ido - i;
            tr2 = cc[i - 1 + (3*k + 2)*ido] + cc[ic - 1 + (3*k + 1)*ido];
            cr2 = cc[i - 1 + 3*k*ido] + taur*tr2;
            ch[i - 1 + k*ido] = cc[i - 1 + 3*k*ido] + tr2;
            ti2 = cc[i + (3*k + 2)*ido]- cc[ic + (3*k + 1)*ido];
            ci2 = cc[i + 3*k*ido] + taur*ti2;
            ch[i + k*ido] = cc[i + 3*k*ido] + ti2;
            cr3 = taui*(cc[i - 1 + (3*k + 2)*ido] - cc[ic - 1 + (3*k + 1)*ido]);
            ci3 = taui*(cc[i + (3*k + 2)*ido] + cc[ic + (3*k + 1)*ido]);
            dr2 = cr2 - ci3;
            dr3 = cr2 + ci3;
            di2 = ci2 + cr3;
            di3 = ci2 - cr3;
            ch[i - 1 + (k + l1)*ido] = wa1[i - 2]*dr2 - wa1[i - 1]*di2;
            ch[i + (k + l1)*ido] = wa1[i - 2]*di2 + wa1[i - 1]*dr2;
            ch[i - 1 + (k + 2*l1)*ido] = wa2[i - 2]*dr3 - wa2[i - 1]*di3;
            ch[i + (k + 2*l1)*ido] = wa2[i - 2]*di3 + wa2[i - 1]*dr3;
        }
    }
} 


static void 
fftpack_radf4(int          ido, 
              int          l1, 
              real   cc[],
              real   ch[],
              real   wa1[],
              real   wa2[],
              real   wa3[])
{
    const real hsqt2 = 0.7071067811865475;
    int i, k, ic;
    real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    
    for (k=0; k<l1; k++) 
    {
        tr1 = cc[(k + l1)*ido] + cc[(k + 3*l1)*ido];
        tr2 = cc[k*ido] + cc[(k + 2*l1)*ido];
        ch[4*k*ido] = tr1 + tr2;
        ch[ido-1 + (4*k + 3)*ido] = tr2 - tr1;
        ch[ido-1 + (4*k + 1)*ido] = cc[k*ido] - cc[(k + 2*l1)*ido];
        ch[(4*k + 2)*ido] = cc[(k + 3*l1)*ido] - cc[(k + l1)*ido];
    }
    if (ido < 2) 
        return;
    if (ido != 2)
    {
        for (k=0; k<l1; k++) 
        {
            for (i=2; i<ido; i += 2)
            {
                ic = ido - i;
                cr2 = wa1[i - 2]*cc[i - 1 + (k + l1)*ido] + wa1[i - 1]*cc[i + (k + l1)*ido];
                ci2 = wa1[i - 2]*cc[i + (k + l1)*ido] - wa1[i - 1]*cc[i - 1 + (k + l1)*ido];
                cr3 = wa2[i - 2]*cc[i - 1 + (k + 2*l1)*ido] + wa2[i - 1]*cc[i + (k + 2*l1)*ido];
                ci3 = wa2[i - 2]*cc[i + (k + 2*l1)*ido] - wa2[i - 1]*cc[i - 1 + (k + 2*l1)*ido];
                cr4 = wa3[i - 2]*cc[i - 1 + (k + 3*l1)*ido] + wa3[i - 1]*cc[i + (k + 3*l1)*ido];
                ci4 = wa3[i - 2]*cc[i + (k + 3*l1)*ido] - wa3[i - 1]*cc[i - 1 + (k + 3*l1)*ido];
                tr1 = cr2 + cr4;
                tr4 = cr4 - cr2;
                ti1 = ci2 + ci4;
                ti4 = ci2 - ci4;
                ti2 = cc[i + k*ido] + ci3;
                ti3 = cc[i + k*ido] - ci3;
                tr2 = cc[i - 1 + k*ido] + cr3;
                tr3 = cc[i - 1 + k*ido] - cr3;
                ch[i - 1 + 4*k*ido] = tr1 + tr2;
                ch[ic - 1 + (4*k + 3)*ido] = tr2 - tr1;
                ch[i + 4*k*ido] = ti1 + ti2;
                ch[ic + (4*k + 3)*ido] = ti1 - ti2;
                ch[i - 1 + (4*k + 2)*ido] = ti4 + tr3;
                ch[ic - 1 + (4*k + 1)*ido] = tr3 - ti4;
                ch[i + (4*k + 2)*ido] = tr4 + ti3;
                ch[ic + (4*k + 1)*ido] = tr4 - ti3;
            }
        }
        if (ido % 2 == 1) 
            return;
    }
    for (k=0; k<l1; k++) 
    {
        ti1 = -hsqt2*(cc[ido-1 + (k + l1)*ido] + cc[ido-1 + (k + 3*l1)*ido]);
        tr1 = hsqt2*(cc[ido-1 + (k + l1)*ido] - cc[ido-1 + (k + 3*l1)*ido]);
        ch[ido-1 + 4*k*ido] = tr1 + cc[ido-1 + k*ido];
        ch[ido-1 + (4*k + 2)*ido] = cc[ido-1 + k*ido] - tr1;
        ch[(4*k + 1)*ido] = ti1 - cc[ido-1 + (k + 2*l1)*ido];
        ch[(4*k + 3)*ido] = ti1 + cc[ido-1 + (k + 2*l1)*ido];
    }
} 


static void
fftpack_radb4(int          ido,
              int          l1,
              real   cc[], 
              real   ch[],
              real   wa1[], 
              real   wa2[], 
              real   wa3[])
{
    const real sqrt2 = 1.414213562373095;
    int i, k, ic;
    real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    for (k = 0; k < l1; k++)
    {
        tr1 = cc[4*k*ido] - cc[ido-1 + (4*k + 3)*ido];
        tr2 = cc[4*k*ido] + cc[ido-1 + (4*k + 3)*ido];
        tr3 = cc[ido-1 + (4*k + 1)*ido] + cc[ido-1 + (4*k + 1)*ido];
        tr4 = cc[(4*k + 2)*ido] + cc[(4*k + 2)*ido];
        ch[k*ido] = tr2 + tr3;
        ch[(k + l1)*ido] = tr1 - tr4;
        ch[(k + 2*l1)*ido] = tr2 - tr3;
        ch[(k + 3*l1)*ido] = tr1 + tr4;
    }
    if (ido < 2) 
        return;
    if (ido != 2) 
    {
        for (k = 0; k < l1; ++k)
        {
            for (i = 2; i < ido; i += 2)
            {
                ic = ido - i;
                ti1 = cc[i + 4*k*ido] + cc[ic + (4*k + 3)*ido];
                ti2 = cc[i + 4*k*ido] - cc[ic + (4*k + 3)*ido];
                ti3 = cc[i + (4*k + 2)*ido] - cc[ic + (4*k + 1)*ido];
                tr4 = cc[i + (4*k + 2)*ido] + cc[ic + (4*k + 1)*ido];
                tr1 = cc[i - 1 + 4*k*ido] - cc[ic - 1 + (4*k + 3)*ido];
                tr2 = cc[i - 1 + 4*k*ido] + cc[ic - 1 + (4*k + 3)*ido];
                ti4 = cc[i - 1 + (4*k + 2)*ido] - cc[ic - 1 + (4*k + 1)*ido];
                tr3 = cc[i - 1 + (4*k + 2)*ido] + cc[ic - 1 + (4*k + 1)*ido];
                ch[i - 1 + k*ido] = tr2 + tr3;
                cr3 = tr2 - tr3;
                ch[i + k*ido] = ti2 + ti3;
                ci3 = ti2 - ti3;
                cr2 = tr1 - tr4;
                cr4 = tr1 + tr4;
                ci2 = ti1 + ti4;
                ci4 = ti1 - ti4;
                ch[i - 1 + (k + l1)*ido] = wa1[i - 2]*cr2 - wa1[i - 1]*ci2;
                ch[i + (k + l1)*ido] = wa1[i - 2]*ci2 + wa1[i - 1]*cr2;
                ch[i - 1 + (k + 2*l1)*ido] = wa2[i - 2]*cr3 - wa2[i - 1]*ci3;
                ch[i + (k + 2*l1)*ido] = wa2[i - 2]*ci3 + wa2[i - 1]*cr3;
                ch[i - 1 + (k + 3*l1)*ido] = wa3[i - 2]*cr4 - wa3[i - 1]*ci4;
                ch[i + (k + 3*l1)*ido] = wa3[i - 2]*ci4 + wa3[i - 1]*cr4;
            }
        }
        if (ido % 2 == 1)
            return;
    }
    for (k = 0; k < l1; k++) 
    {
        ti1 = cc[(4*k + 1)*ido] + cc[(4*k + 3)*ido];
        ti2 = cc[(4*k + 3)*ido] - cc[(4*k + 1)*ido];
        tr1 = cc[ido-1 + 4*k*ido] - cc[ido-1 + (4*k + 2)*ido];
        tr2 = cc[ido-1 + 4*k*ido] + cc[ido-1 + (4*k + 2)*ido];
        ch[ido-1 + k*ido] = tr2 + tr2;
        ch[ido-1 + (k + l1)*ido] = sqrt2*(tr1 - ti1);
        ch[ido-1 + (k + 2*l1)*ido] = ti2 + ti2;
        ch[ido-1 + (k + 3*l1)*ido] = -sqrt2*(tr1 + ti1);
    }
}


static void 
fftpack_radf5(int          ido,
              int          l1, 
              real   cc[], 
              real   ch[],
              real   wa1[],
              real   wa2[], 
              real   wa3[], 
              real   wa4[])
{
    const real tr11 = 0.309016994374947;
    const real ti11 = 0.951056516295154;
    const real tr12 = -0.809016994374947;
    const real ti12 = 0.587785252292473;
    int i, k, ic;
    real ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, dr4, dr5,
        cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
 
    for (k = 0; k < l1; k++) 
    {
        cr2 = cc[(k + 4*l1)*ido] + cc[(k + l1)*ido];
        ci5 = cc[(k + 4*l1)*ido] - cc[(k + l1)*ido];
        cr3 = cc[(k + 3*l1)*ido] + cc[(k + 2*l1)*ido];
        ci4 = cc[(k + 3*l1)*ido] - cc[(k + 2*l1)*ido];
        ch[5*k*ido] = cc[k*ido] + cr2 + cr3;
        ch[ido-1 + (5*k + 1)*ido] = cc[k*ido] + tr11*cr2 + tr12*cr3;
        ch[(5*k + 2)*ido] = ti11*ci5 + ti12*ci4;
        ch[ido-1 + (5*k + 3)*ido] = cc[k*ido] + tr12*cr2 + tr11*cr3;
        ch[(5*k + 4)*ido] = ti12*ci5 - ti11*ci4;
    }
    if (ido == 1) 
        return;
    for (k = 0; k < l1; ++k) 
    {
        for (i = 2; i < ido; i += 2)
        {
            ic = ido - i;
            dr2 = wa1[i - 2]*cc[i - 1 + (k + l1)*ido] + wa1[i - 1]*cc[i + (k + l1)*ido];
            di2 = wa1[i - 2]*cc[i + (k + l1)*ido] - wa1[i - 1]*cc[i - 1 + (k + l1)*ido];
            dr3 = wa2[i - 2]*cc[i - 1 + (k + 2*l1)*ido] + wa2[i - 1]*cc[i + (k + 2*l1)*ido];
            di3 = wa2[i - 2]*cc[i + (k + 2*l1)*ido] - wa2[i - 1]*cc[i - 1 + (k + 2*l1)*ido];
            dr4 = wa3[i - 2]*cc[i - 1 + (k + 3*l1)*ido] + wa3[i - 1]*cc[i + (k + 3*l1)*ido];
            di4 = wa3[i - 2]*cc[i + (k + 3*l1)*ido] - wa3[i - 1]*cc[i - 1 + (k + 3*l1)*ido];
            dr5 = wa4[i - 2]*cc[i - 1 + (k + 4*l1)*ido] + wa4[i - 1]*cc[i + (k + 4*l1)*ido];
            di5 = wa4[i - 2]*cc[i + (k + 4*l1)*ido] - wa4[i - 1]*cc[i - 1 + (k + 4*l1)*ido];
            cr2 = dr2 + dr5;
            ci5 = dr5 - dr2;
            cr5 = di2 - di5;
            ci2 = di2 + di5;
            cr3 = dr3 + dr4;
            ci4 = dr4 - dr3;
            cr4 = di3 - di4;
            ci3 = di3 + di4;
            ch[i - 1 + 5*k*ido] = cc[i - 1 + k*ido] + cr2 + cr3;
            ch[i + 5*k*ido] = cc[i + k*ido] + ci2 + ci3;
            tr2 = cc[i - 1 + k*ido] + tr11*cr2 + tr12*cr3;
            ti2 = cc[i + k*ido] + tr11*ci2 + tr12*ci3;
            tr3 = cc[i - 1 + k*ido] + tr12*cr2 + tr11*cr3;
            ti3 = cc[i + k*ido] + tr12*ci2 + tr11*ci3;
            tr5 = ti11*cr5 + ti12*cr4;
            ti5 = ti11*ci5 + ti12*ci4;
            tr4 = ti12*cr5 - ti11*cr4;
            ti4 = ti12*ci5 - ti11*ci4;
            ch[i - 1 + (5*k + 2)*ido] = tr2 + tr5;
            ch[ic - 1 + (5*k + 1)*ido] = tr2 - tr5;
            ch[i + (5*k + 2)*ido] = ti2 + ti5;
            ch[ic + (5*k + 1)*ido] = ti5 - ti2;
            ch[i - 1 + (5*k + 4)*ido] = tr3 + tr4;
            ch[ic - 1 + (5*k + 3)*ido] = tr3 - tr4;
            ch[i + (5*k + 4)*ido] = ti3 + ti4;
            ch[ic + (5*k + 3)*ido] = ti4 - ti3;
        }
    }
} 


static void
fftpack_radb5(int          ido,
              int          l1, 
              real   cc[], 
              real   ch[],
              real   wa1[],
              real   wa2[], 
              real   wa3[],
              real   wa4[])
{
    const real tr11 = 0.309016994374947;
    const real ti11 = 0.951056516295154;
    const real tr12 = -0.809016994374947;
    const real ti12 = 0.587785252292473;
    
    int i, k, ic;
    real ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
        ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    
    for (k = 0; k < l1; k++)
    {
        ti5 = 2*cc[(5*k + 2)*ido];
        ti4 = 2*cc[(5*k + 4)*ido];
        tr2 = 2*cc[ido-1 + (5*k + 1)*ido];
        tr3 = 2*cc[ido-1 + (5*k + 3)*ido];
        ch[k*ido] = cc[5*k*ido] + tr2 + tr3;
        cr2 = cc[5*k*ido] + tr11*tr2 + tr12*tr3;
        cr3 = cc[5*k*ido] + tr12*tr2 + tr11*tr3;
        ci5 = ti11*ti5 + ti12*ti4;
        ci4 = ti12*ti5 - ti11*ti4;
        ch[(k + l1)*ido] = cr2 - ci5;
        ch[(k + 2*l1)*ido] = cr3 - ci4;
        ch[(k + 3*l1)*ido] = cr3 + ci4;
        ch[(k + 4*l1)*ido] = cr2 + ci5;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; ++k)
    {
        for (i = 2; i < ido; i += 2)
        {
            ic = ido - i;
            ti5 = cc[i + (5*k + 2)*ido] + cc[ic + (5*k + 1)*ido];
            ti2 = cc[i + (5*k + 2)*ido] - cc[ic + (5*k + 1)*ido];
            ti4 = cc[i + (5*k + 4)*ido] + cc[ic + (5*k + 3)*ido];
            ti3 = cc[i + (5*k + 4)*ido] - cc[ic + (5*k + 3)*ido];
            tr5 = cc[i - 1 + (5*k + 2)*ido] - cc[ic - 1 + (5*k + 1)*ido];
            tr2 = cc[i - 1 + (5*k + 2)*ido] + cc[ic - 1 + (5*k + 1)*ido];
            tr4 = cc[i - 1 + (5*k + 4)*ido] - cc[ic - 1 + (5*k + 3)*ido];
            tr3 = cc[i - 1 + (5*k + 4)*ido] + cc[ic - 1 + (5*k + 3)*ido];
            ch[i - 1 + k*ido] = cc[i - 1 + 5*k*ido] + tr2 + tr3;
            ch[i + k*ido] = cc[i + 5*k*ido] + ti2 + ti3;
            cr2 = cc[i - 1 + 5*k*ido] + tr11*tr2 + tr12*tr3;
            ci2 = cc[i + 5*k*ido] + tr11*ti2 + tr12*ti3;
            cr3 = cc[i - 1 + 5*k*ido] + tr12*tr2 + tr11*tr3;
            ci3 = cc[i + 5*k*ido] + tr12*ti2 + tr11*ti3;
            cr5 = ti11*tr5 + ti12*tr4;
            ci5 = ti11*ti5 + ti12*ti4;
            cr4 = ti12*tr5 - ti11*tr4;
            ci4 = ti12*ti5 - ti11*ti4;
            dr3 = cr3 - ci4;
            dr4 = cr3 + ci4;
            di3 = ci3 + cr4;
            di4 = ci3 - cr4;
            dr5 = cr2 + ci5;
            dr2 = cr2 - ci5;
            di5 = ci2 - cr5;
            di2 = ci2 + cr5;
            ch[i - 1 + (k + l1)*ido] = wa1[i - 2]*dr2 - wa1[i - 1]*di2;
            ch[i + (k + l1)*ido] = wa1[i - 2]*di2 + wa1[i - 1]*dr2;
            ch[i - 1 + (k + 2*l1)*ido] = wa2[i - 2]*dr3 - wa2[i - 1]*di3;
            ch[i + (k + 2*l1)*ido] = wa2[i - 2]*di3 + wa2[i - 1]*dr3;
            ch[i - 1 + (k + 3*l1)*ido] = wa3[i - 2]*dr4 - wa3[i - 1]*di4;
            ch[i + (k + 3*l1)*ido] = wa3[i - 2]*di4 + wa3[i - 1]*dr4;
            ch[i - 1 + (k + 4*l1)*ido] = wa4[i - 2]*dr5 - wa4[i - 1]*di5;
            ch[i + (k + 4*l1)*ido] = wa4[i - 2]*di5 + wa4[i - 1]*dr5;
        }
    }
} 


static void 
fftpack_radfg(int          ido, 
              int          ip,
              int          l1, 
              int          idl1,
              real   cc[], 
              real   ch[],
              real   wa[])
{
    const real twopi = 6.28318530717959;
    int idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is, nbd;
    real dc2, ai1, ai2, ar1, ar2, ds2, dcp, arg, dsp, ar1h, ar2h;
    arg = twopi / ip;
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (ip + 1) / 2;
    nbd = (ido - 1) / 2;
    if (ido != 1)
    {
        for (ik=0; ik<idl1; ik++) ch[ik] = cc[ik];
        for (j=1; j<ip; j++)
            for (k=0; k<l1; k++)
                ch[(k + j*l1)*ido] = cc[(k + j*l1)*ido];
        if (nbd <= l1) 
        {
            is = -ido;
            for (j=1; j<ip; j++) 
            {
                is += ido;
                idij = is-1;
                for (i=2; i<ido; i+=2) 
                {
                    idij += 2;
                    for (k=0; k<l1; k++) 
                    {
                        ch[i - 1 + (k + j*l1)*ido] =
                        wa[idij - 1]*cc[i - 1 + (k + j*l1)*ido] + wa[idij]*cc[i + (k + j*l1)*ido];
                        ch[i + (k + j*l1)*ido] =
                            wa[idij - 1]*cc[i + (k + j*l1)*ido] - wa[idij]*cc[i - 1 + (k + j*l1)*ido];
                    }
                }
            }
        }
        else 
        {
            is = -ido;
            for (j=1; j<ip; j++)
            {
                is += ido;
                for (k=0; k<l1; k++)
                {
                    idij = is-1;
                    for (i=2; i<ido; i+=2) 
                    {
                        idij += 2;
                        ch[i - 1 + (k + j*l1)*ido] =
                            wa[idij - 1]*cc[i - 1 + (k + j*l1)*ido] + wa[idij]*cc[i + (k + j*l1)*ido];
                        ch[i + (k + j*l1)*ido] =
                            wa[idij - 1]*cc[i + (k + j*l1)*ido] - wa[idij]*cc[i - 1 + (k + j*l1)*ido];
                    }
                }
            }
        }
        if (nbd >= l1) 
        {
            for (j=1; j<ipph; j++) 
            {
                jc = ip - j;
                for (k=0; k<l1; k++)
                {
                    for (i=2; i<ido; i+=2) 
                    {
                        cc[i - 1 + (k + j*l1)*ido] = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
                        cc[i - 1 + (k + jc*l1)*ido] = ch[i + (k + j*l1)*ido] - ch[i + (k + jc*l1)*ido];
                        cc[i + (k + j*l1)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
                        cc[i + (k + jc*l1)*ido] = ch[i - 1 + (k + jc*l1)*ido] - ch[i - 1 + (k + j*l1)*ido];
                    }
                }
            }
        }
        else 
        {
            for (j=1; j<ipph; j++) 
            {
                jc = ip - j;
                for (i=2; i<ido; i+=2)
                {
                    for (k=0; k<l1; k++)
                    {
                        cc[i - 1 + (k + j*l1)*ido] =
                        ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
                        cc[i - 1 + (k + jc*l1)*ido] = ch[i + (k + j*l1)*ido] - ch[i + (k + jc*l1)*ido];
                        cc[i + (k + j*l1)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
                        cc[i + (k + jc*l1)*ido] = ch[i - 1 + (k + jc*l1)*ido] - ch[i - 1 + (k + j*l1)*ido];
                    }
                }
            }
        }
    }
    else
    {
        for (ik=0; ik<idl1; ik++) 
            cc[ik] = ch[ik];
    }
    for (j=1; j<ipph; j++)
    {
        jc = ip - j;
        for (k=0; k<l1; k++) 
        {
            cc[(k + j*l1)*ido] = ch[(k + j*l1)*ido] + ch[(k + jc*l1)*ido];
            cc[(k + jc*l1)*ido] = ch[(k + jc*l1)*ido] - ch[(k + j*l1)*ido];
        }
    }
    
    ar1 = 1;
    ai1 = 0;
    for (l=1; l<ipph; l++) 
    {
        lc = ip - l;
        ar1h = dcp*ar1 - dsp*ai1;
        ai1 = dcp*ai1 + dsp*ar1;
        ar1 = ar1h;
        for (ik=0; ik<idl1; ik++) 
        {
            ch[ik + l*idl1] = cc[ik] + ar1*cc[ik + idl1];
            ch[ik + lc*idl1] = ai1*cc[ik + (ip-1)*idl1];
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        for (j=2; j<ipph; j++)
        {
            jc = ip - j;
            ar2h = dc2*ar2 - ds2*ai2;
            ai2 = dc2*ai2 + ds2*ar2;
            ar2 = ar2h;
            for (ik=0; ik<idl1; ik++) 
            {
                ch[ik + l*idl1] += ar2*cc[ik + j*idl1];
                ch[ik + lc*idl1] += ai2*cc[ik + jc*idl1];
            }
        }
    }
    for (j=1; j<ipph; j++)
        for (ik=0; ik<idl1; ik++)
            ch[ik] += cc[ik + j*idl1];
    
    if (ido >= l1) 
    {
        for (k=0; k<l1; k++) 
        {
            for (i=0; i<ido; i++) 
            {
                cc[i + k*ip*ido] = ch[i + k*ido];
            }
        }
    } 
    else
    {
        for (i=0; i<ido; i++)
        {
            for (k=0; k<l1; k++)
            {
                cc[i + k*ip*ido] = ch[i + k*ido];
            }
        }
    }
    for (j=1; j<ipph; j++)
    {
        jc = ip - j;
        j2 = 2*j;
        for (k=0; k<l1; k++)
        {
            cc[ido-1 + (j2 - 1 + k*ip)*ido] = ch[(k + j*l1)*ido];
            cc[(j2 + k*ip)*ido] = ch[(k + jc*l1)*ido];
        }
    }
    if (ido == 1) return;
    if (nbd >= l1)
    {
        for (j=1; j<ipph; j++) 
        {
            jc = ip - j;
            j2 = 2*j;
            for (k=0; k<l1; k++)
            {
                for (i=2; i<ido; i+=2)
                {
                    ic = ido - i;
                    cc[i - 1 + (j2 + k*ip)*ido] = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
                    cc[ic - 1 + (j2 - 1 + k*ip)*ido] = ch[i - 1 + (k + j*l1)*ido] - ch[i - 1 + (k + jc*l1)*ido];
                    cc[i + (j2 + k*ip)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
                    cc[ic + (j2 - 1 + k*ip)*ido] = ch[i + (k + jc*l1)*ido] - ch[i + (k + j*l1)*ido];
                }
            }
        }
    }
    else
    {
        for (j=1; j<ipph; j++)
        {
            jc = ip - j;
            j2 = 2*j;
            for (i=2; i<ido; i+=2) 
            {
                ic = ido - i;
                for (k=0; k<l1; k++)
                {
                    cc[i - 1 + (j2 + k*ip)*ido] = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
                    cc[ic - 1 + (j2 - 1 + k*ip)*ido] = ch[i - 1 + (k + j*l1)*ido] - ch[i - 1 + (k + jc*l1)*ido];
                    cc[i + (j2 + k*ip)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
                    cc[ic + (j2 - 1 + k*ip)*ido] = ch[i + (k + jc*l1)*ido] - ch[i + (k + j*l1)*ido];
                }
            }
        }
    }
} 


static void 
fftpack_radbg(int          ido,
              int          ip, 
              int          l1, 
              int          idl1,
              real   cc[],
              real   ch[],
              real   wa[])
{
    const real twopi = 6.28318530717959;
    int idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is;
    real dc2, ai1, ai2, ar1, ar2, ds2;
    int nbd;
    real dcp, arg, dsp, ar1h, ar2h;
    arg = twopi / ip;
    dcp = cos(arg);
    dsp = sin(arg);
    nbd = (ido - 1) / 2;
    ipph = (ip + 1) / 2;
    
    if (ido >= l1) 
    {
        for (k=0; k<l1; k++) 
        {
            for (i=0; i<ido; i++)
            {
                ch[i + k*ido] = cc[i + k*ip*ido];
            }
        }
    }
    else 
    {
        for (i=0; i<ido; i++) 
        {
            for (k=0; k<l1; k++)
            {
                ch[i + k*ido] = cc[i + k*ip*ido];
            }
        }
    }
    for (j=1; j<ipph; j++)
    {
        jc = ip - j;
        j2 = 2*j;
        for (k=0; k<l1; k++) 
        {
            ch[(k + j*l1)*ido] = cc[ido-1 + (j2 - 1 + k*ip)*ido] + cc[ido-1 + (j2 - 1 + k*ip)*ido];
            ch[(k + jc*l1)*ido] = cc[(j2 + k*ip)*ido] + cc[(j2 + k*ip)*ido];
        }
    }
    
    if (ido != 1) 
    {
        if (nbd >= l1)
        {
            for (j=1; j<ipph; j++) 
            {
                jc = ip - j;
                for (k=0; k<l1; k++)
                {
                    for (i=2; i<ido; i+=2) 
                    {
                        ic = ido - i;
                        ch[i - 1 + (k + j*l1)*ido] = cc[i - 1 + (2*j + k*ip)*ido] +  cc[ic - 1 + (2*j - 1 + k*ip)*ido];
                        ch[i - 1 + (k + jc*l1)*ido] = cc[i - 1 + (2*j + k*ip)*ido] - cc[ic - 1 + (2*j - 1 + k*ip)*ido];
                        ch[i + (k + j*l1)*ido] = cc[i + (2*j + k*ip)*ido] - cc[ic + (2*j - 1 + k*ip)*ido];
                        ch[i + (k + jc*l1)*ido] = cc[i + (2*j + k*ip)*ido] + cc[ic + (2*j - 1 + k*ip)*ido];
                    }
                }
            }
        }
        else 
        {
            for (j=1; j<ipph; j++) 
            {
                jc = ip - j;
                for (i=2; i<ido; i+=2) 
                {
                    ic = ido - i;
                    for (k=0; k<l1; k++) 
                    {
                        ch[i - 1 + (k + j*l1)*ido] = cc[i - 1 + (2*j + k*ip)*ido] + cc[ic - 1 + (2*j - 1 + k*ip)*ido];
                        ch[i - 1 + (k + jc*l1)*ido] = cc[i - 1 + (2*j + k*ip)*ido] - cc[ic - 1 + (2*j - 1 + k*ip)*ido];
                        ch[i + (k + j*l1)*ido] = cc[i + (2*j + k*ip)*ido] - cc[ic + (2*j - 1 + k*ip)*ido];
                        ch[i + (k + jc*l1)*ido] = cc[i + (2*j + k*ip)*ido] + cc[ic + (2*j - 1 + k*ip)*ido];
                    }
                }
            }
        }
    }
    
    ar1 = 1;
    ai1 = 0;
    for (l=1; l<ipph; l++)
    {
        lc = ip - l;
        ar1h = dcp*ar1 - dsp*ai1;
        ai1 = dcp*ai1 + dsp*ar1;
        ar1 = ar1h;
        for (ik=0; ik<idl1; ik++) 
        {
            cc[ik + l*idl1] = ch[ik] + ar1*ch[ik + idl1];
            cc[ik + lc*idl1] = ai1*ch[ik + (ip-1)*idl1];
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        for (j=2; j<ipph; j++)
        {
            jc = ip - j;
            ar2h = dc2*ar2 - ds2*ai2;
            ai2 = dc2*ai2 + ds2*ar2;
            ar2 = ar2h;
            for (ik=0; ik<idl1; ik++) 
            {
                cc[ik + l*idl1] += ar2*ch[ik + j*idl1];
                cc[ik + lc*idl1] += ai2*ch[ik + jc*idl1];
            }
        }
    }
    for (j=1; j<ipph; j++) 
    {
        for (ik=0; ik<idl1; ik++) 
        {
            ch[ik] += ch[ik + j*idl1];
        }
    }
    for (j=1; j<ipph; j++) 
    {
        jc = ip - j;
        for (k=0; k<l1; k++) 
        {
            ch[(k + j*l1)*ido] = cc[(k + j*l1)*ido] - cc[(k + jc*l1)*ido];
            ch[(k + jc*l1)*ido] = cc[(k + j*l1)*ido] + cc[(k + jc*l1)*ido];
        }
    }
    
    if (ido == 1) return;
    if (nbd >= l1) 
    {
        for (j=1; j<ipph; j++) 
        {
            jc = ip - j;
            for (k=0; k<l1; k++) 
            {
                for (i=2; i<ido; i+=2)
                {
                    ch[i - 1 + (k + j*l1)*ido] = cc[i - 1 + (k + j*l1)*ido] - cc[i + (k + jc*l1)*ido];
                    ch[i - 1 + (k + jc*l1)*ido] = cc[i - 1 + (k + j*l1)*ido] + cc[i + (k + jc*l1)*ido];
                    ch[i + (k + j*l1)*ido] = cc[i + (k + j*l1)*ido] + cc[i - 1 + (k + jc*l1)*ido];
                    ch[i + (k + jc*l1)*ido] = cc[i + (k + j*l1)*ido] - cc[i - 1 + (k + jc*l1)*ido];
                }
            }
        }
    }
    else
    {
        for (j=1; j<ipph; j++) 
        {
            jc = ip - j;
            for (i=2; i<ido; i+=2)
            {
                for (k=0; k<l1; k++)
                {
                    ch[i - 1 + (k + j*l1)*ido] = cc[i - 1 + (k + j*l1)*ido] - cc[i + (k + jc*l1)*ido];
                    ch[i - 1 + (k + jc*l1)*ido] = cc[i - 1 + (k + j *l1)*ido] + cc[i + (k + jc*l1)*ido];
                    ch[i + (k + j*l1)*ido] = cc[i + (k + j*l1)*ido] + cc[i - 1 + (k + jc*l1)*ido];
                    ch[i + (k + jc*l1)*ido] = cc[i + (k + j*l1)*ido] - cc[i - 1 + (k + jc*l1)*ido];
                }
            }
        }
    }
    for (ik=0; ik<idl1; ik++)
    {
        cc[ik] = ch[ik];
    }
    for (j=1; j<ip; j++)
        for (k=0; k<l1; k++)
            cc[(k + j*l1)*ido] = ch[(k + j*l1)*ido];

    if (nbd <= l1) 
    {
        is = -ido;
        for (j=1; j<ip; j++)
        {
            is += ido;
            idij = is-1;
            for (i=2; i<ido; i+=2)
            {
                idij += 2;
                for (k=0; k<l1; k++) 
                {
                    cc[i - 1 + (k + j*l1)*ido] = wa[idij - 1]*ch[i - 1 + (k + j*l1)*ido] - wa[idij]*ch[i + (k + j*l1)*ido];
                    cc[i + (k + j*l1)*ido] = wa[idij - 1]*ch[i + (k + j*l1)*ido] + wa[idij]*ch[i - 1 + (k + j*l1)*ido];
                }
            }
        }
    }
    else
    {
        is = -ido;
        for (j=1; j<ip; j++) 
        {
            is += ido;
            for (k=0; k<l1; k++)
            {
                idij = is;
                for (i=2; i<ido; i+=2)
                {
                    idij += 2;
                    cc[i - 1 + (k + j*l1)*ido] = wa[idij-1]*ch[i - 1 + (k + j*l1)*ido] - wa[idij]*ch[i + (k + j*l1)*ido];
                    cc[i + (k + j*l1)*ido] = wa[idij-1]*ch[i + (k + j*l1)*ido] + wa[idij]*ch[i - 1 + (k + j*l1)*ido];
                }
            }
        }
    }
} 



static void 
fftpack_cfftf1(int          n,
               real   c[],
               real   ch[],
               real   wa[],
               int          ifac[15],
               int          isign)
{
    int idot, i;
    int k1, l1, l2;
    int na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;
    real *cinput, *coutput;
    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 0;
    
    for (k1=2; k1<=nf+1; k1++) 
    {
        ip = ifac[k1];
        l2 = ip*l1;
        ido = n / l2;
        idot = ido + ido;
        idl1 = idot*l1;
        if (na) 
        {
            cinput = ch;
            coutput = c;
        }
        else 
        {
            cinput = c;
            coutput = ch;
        }
        switch (ip) 
        {
            case 4:
                ix2 = iw + idot;
                ix3 = ix2 + idot;
                fftpack_passf4(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], isign);
                na = !na;
                break;
            case 2:
                fftpack_passf2(idot, l1, cinput, coutput, &wa[iw], isign);
                na = !na;
                break;
            case 3:
                ix2 = iw + idot;
                fftpack_passf3(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], isign);
                na = !na;
                break;
            case 5:
                ix2 = iw + idot;
                ix3 = ix2 + idot;
                ix4 = ix3 + idot;
                fftpack_passf5(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4], isign);
                na = !na;
                break;
            default:
                fftpack_passf(&nac, idot, ip, l1, idl1, cinput, coutput, &wa[iw], isign);
                if (nac != 0) na = !na;
        }
        l1 = l2;
        iw += (ip - 1)*idot;
    }
    if (na == 0) 
        return;
    for (i=0; i<2*n; i++) 
        c[i] = ch[i];
}


void
fftpack_cfftf(int          n, 
              real   c[], 
              real   wsave[])
{
    int iw1, iw2;
    
    if (n == 1)
        return;
    iw1 = 2*n;
    iw2 = iw1 + 2*n;
    fftpack_cfftf1(n, c, wsave, wsave+iw1, (int*)(wsave+iw2), -1);
} 


void 
fftpack_cfftb(int          n, 
              real   c[], 
              real   wsave[])
{
    int iw1, iw2;
    
    if (n == 1)
        return;
    iw1 = 2*n;
    iw2 = iw1 + 2*n;
    fftpack_cfftf1(n, c, wsave, wsave+iw1, (int*)(wsave+iw2), +1);
} 


static void 
fftpack_factorize(int    n,
                  int    ifac[15])
{
    static const int ntryh[4] = { 3,4,2,5 };
    int ntry=3, i, j=0, ib, nf=0, nl=n, nq, nr;

startloop:
    if (j < 4)
        ntry = ntryh[j];
    else
        ntry+= 2;
    j++;
    do 
    {
        nq = nl / ntry;
        nr = nl - ntry*nq;
        if (nr != 0) goto startloop;
        nf++;
        ifac[nf + 1] = ntry;
        nl = nq;
        if (ntry == 2 && nf != 1) 
        {
            for (i=2; i<=nf; i++) 
            {
                ib = nf - i + 2;
                ifac[ib + 1] = ifac[ib];
            }
            ifac[2] = 2;
        }
    } 
    while (nl != 1);
    ifac[0] = n;
    ifac[1] = nf;
}


static void
fftpack_cffti1(int          n, 
               real   wa[], 
               int          ifac[15])
{
    const real twopi = 6.28318530717959;
    real arg, argh, argld, fi;
    int idot, i, j;
    int i1, k1, l1, l2;
    int ld, ii, nf, ip;
    int ido, ipm;

    fftpack_factorize(n,ifac);
    nf = ifac[1];
    argh = twopi/(real)n;
    i = 1;
    l1 = 1;
    for (k1=1; k1<=nf; k1++) 
    {
        ip = ifac[k1+1];
        ld = 0;
        l2 = l1*ip;
        ido = n / l2;
        idot = ido + ido + 2;
        ipm = ip - 1;
        for (j=1; j<=ipm; j++)
        {
            i1 = i;
            wa[i-1] = 1;
            wa[i] = 0;
            ld += l1;
            fi = 0;
            argld = ld*argh;
            for (ii=4; ii<=idot; ii+=2) 
            {
                i+= 2;
                fi+= 1;
                arg = fi*argld;
                wa[i-1] = cos(arg);
                wa[i] = sin(arg);
            }
            if (ip > 5) 
            {
                wa[i1-1] = wa[i-1];
                wa[i1] = wa[i];
            }
        }
        l1 = l2;
    }
} 




static void 
fftpack_rfftf1(int n, 
               real   c[],
               real   ch[], 
               real   wa[], 
               int          ifac[15])
{
    int i;
    int k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;
    real *cinput, *coutput;
    nf = ifac[1];
    na = 1;
    l2 = n;
    iw = n-1;
    for (k1 = 1; k1 <= nf; ++k1) 
    {
        kh = nf - k1;
        ip = ifac[kh + 2];
        l1 = l2 / ip;
        ido = n / l2;
        idl1 = ido*l1;
        iw -= (ip - 1)*ido;
        na = !na;
        if (na) 
        {
            cinput = ch;
            coutput = c;
        }
        else 
        {
            cinput = c;
            coutput = ch;
        }
      switch (ip) 
      {
          case 4:
              ix2 = iw + ido;
              ix3 = ix2 + ido;
              fftpack_radf4(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3]);
              break;
          case 2:
              fftpack_radf2(ido, l1, cinput, coutput, &wa[iw]);
              break;
          case 3:
              ix2 = iw + ido;
              fftpack_radf3(ido, l1, cinput, coutput, &wa[iw], &wa[ix2]);
              break;
          case 5:
              ix2 = iw + ido;
              ix3 = ix2 + ido;
              ix4 = ix3 + ido;
              fftpack_radf5(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
              break;
          default:
              if (ido == 1)
                  na = !na;
              if (na == 0)
              {
                  fftpack_radfg(ido, ip, l1, idl1, c, ch, &wa[iw]);
                  na = 1;
              }
                  else 
                  {
                      fftpack_radfg(ido, ip, l1, idl1, ch, c, &wa[iw]);
                      na = 0;
                  }
      }
        l2 = l1;
    }
    if (na == 1)
        return;
    for (i = 0; i < n; i++)
        c[i] = ch[i];
} 


static void 
fftpack_rfftb1(int          n, 
               real   c[],
               real   ch[], 
               real   wa[], 
               int          ifac[15])
{
    int i;
    int k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, idl1;
    real *cinput, *coutput;
    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 0;
    
    for (k1=1; k1<=nf; k1++) 
    {
        ip = ifac[k1 + 1];
        l2 = ip*l1;
        ido = n / l2;
        idl1 = ido*l1;
        if (na) 
        {
            cinput = ch;
            coutput = c;
        }
        else 
        {
            cinput = c;
            coutput = ch;
        }
        switch (ip) 
        {
            case 4:
                ix2 = iw + ido;
                ix3 = ix2 + ido;
                fftpack_radb4(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3]);
                na = !na;
                break;
            case 2:
                fftpack_radb2(ido, l1, cinput, coutput, &wa[iw]);
                na = !na;
                break;
            case 3:
                ix2 = iw + ido;
                fftpack_radb3(ido, l1, cinput, coutput, &wa[iw], &wa[ix2]);
                na = !na;
                break;
            case 5:
                ix2 = iw + ido;
                ix3 = ix2 + ido;
                ix4 = ix3 + ido;
                fftpack_radb5(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
                na = !na;
                break;
            default:
                fftpack_radbg(ido, ip, l1, idl1, cinput, coutput, &wa[iw]);
                if (ido == 1) na = !na;
        }
        l1 = l2;
        iw += (ip - 1)*ido;
    }
    if (na == 0) 
        return;
    for (i=0; i<n; i++) 
        c[i] = ch[i];
} 




static void
fftpack_rffti1(int          n, 
               real         wa[], 
               int          ifac[15])
{
    const real twopi = 6.28318530717959;
    real arg, argh, argld, fi;
    int i, j;
    int k1, l1, l2;
    int ld, ii, nf, ip, is;
    int ido, ipm, nfm1;
    fftpack_factorize(n,ifac);
    nf = ifac[1];
    argh = twopi / n;
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) return;
    for (k1 = 1; k1 <= nfm1; k1++) 
    {
        ip = ifac[k1 + 1];
        ld = 0;
        l2 = l1*ip;
        ido = n / l2;
        ipm = ip - 1;
        for (j = 1; j <= ipm; ++j) 
        {
            ld += l1;
            i = is;
            argld = (real) ld*argh;
            fi = 0;
            for (ii = 3; ii <= ido; ii += 2) 
            {
                i += 2;
                fi += 1;
                arg = fi*argld;
                wa[i - 2] = cos(arg);
                wa[i - 1] = sin(arg);
            }
            is += ido;
        }
        l1 = l2;
    }
} 




/* End of fftpack - begin GROMACS code */


int
gmx_fft_init_1d(gmx_fft_t *        pfft,
                int                nx,
                int                flags)
{
    gmx_fft_t    fft;
   
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    if( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
        
    fft->next = NULL;
    fft->n    = nx;
    
    /* Need 4*n storage for 1D complex FFT */
    if( (fft->work = (real *)malloc(sizeof(real)*(4*nx))) == NULL) 
    {
        free(fft);
        return ENOMEM;
    }

    if(fft->n>1)
        fftpack_cffti1(nx,fft->work,fft->ifac);
    
    *pfft = fft;
    return 0;
};



int
gmx_fft_init_1d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                flags)
{
    gmx_fft_t    fft;
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    

    fft->next = NULL;
    fft->n    = nx;

    /* Need 2*n storage for 1D real FFT */
    if((fft->work = (real *)malloc(sizeof(real)*(2*nx)))==NULL) 
    {
        free(fft);
        return ENOMEM;
    }  
    
    if(fft->n>1)
        fftpack_rffti1(nx,fft->work,fft->ifac);
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_2d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                int                flags)
{
    gmx_fft_t     fft;
    int           rc;
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */
    if( (rc = gmx_fft_init_1d(&fft,nx,flags)) != 0)
    {
        return rc;
    }    
    
    /* Create Y transform as a link from X */
    if( (rc=gmx_fft_init_1d(&(fft->next),ny,flags)) != 0)
    {
        free(fft);
        return rc;
    }
    
    *pfft = fft;
    return 0;
};


int
gmx_fft_init_2d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     int                flags)
{
    gmx_fft_t     fft;
    int           nyc = (ny/2 + 1);
    int           rc;
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    /* Create the X transform */
    if( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    fft->n    = nx;
    
    /* Need 4*nx storage for 1D complex FFT, and another
     * 2*nx*nyc elements for complex-to-real storage in our high-level routine.
     */
    if( (fft->work = (real *)malloc(sizeof(real)*(4*nx+2*nx*nyc))) == NULL) 
    {
        free(fft);
        return ENOMEM;
    }
    fftpack_cffti1(nx,fft->work,fft->ifac);
    
    /* Create real Y transform as a link from X */
    if( (rc=gmx_fft_init_1d_real(&(fft->next),ny,flags)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;
    return 0;
}


int
gmx_fft_init_3d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                int                nz,
                int                flags)
{
    gmx_fft_t     fft;
    int           rc;
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */

    if( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    

    fft->n    = nx;
    
    /* Need 4*nx storage for 1D complex FFT, and another
     * 2*nz elements for gmx_fft_transpose_2d_nelem() storage.
     */
    if( (fft->work = (real *)malloc(sizeof(real)*(4*nx+2*nz))) == NULL) 
    {
        free(fft);
        return ENOMEM;
    }
    
    fftpack_cffti1(nx,fft->work,fft->ifac);

    
    /* Create 2D Y/Z transforms as a link from X */
    if( (rc=gmx_fft_init_2d(&(fft->next),ny,nz,flags)) != 0)
    {
        free(fft);
        return rc;
    }
    
    *pfft = fft;
    return 0;
};


int
gmx_fft_init_3d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     int                nz,
                     int                flags)
{
    gmx_fft_t     fft;
    int           nzc = (nz/2 + 1);
    int           rc;
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;
        
    /* Create the X transform */
    if( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    fft->n    = nx;

    /* Need 4*nx storage for 1D complex FFT, another
     * 2*nx*ny*nzc elements to copy the entire 3D matrix when
     * doing out-of-place complex-to-real FFTs, and finally
     * 2*nzc elements for transpose work space.
     */
    if( (fft->work = (real *)malloc(sizeof(real)*(4*nx+2*nx*ny*nzc+2*nzc))) == NULL) 
    {
        free(fft);
        return ENOMEM;
    }
    fftpack_cffti1(nx,fft->work,fft->ifac);
    
    /* Create 2D real Y/Z transform as a link from X */
    if( (rc=gmx_fft_init_2d_real(&(fft->next),ny,nz,flags)) != 0)
    {
        free(fft);
        return rc;
    }
    
    *pfft = fft;
    return 0;
}


int 
gmx_fft_1d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int             i,n;
    real *    p1;
    real *    p2;

    n=fft->n;

    if(n==1)
    {
        p1 = (real *)in_data;
        p2 = (real *)out_data;        
        p2[0] = p1[0];
        p2[1] = p1[1];
    }
    
    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     */
    if( in_data != out_data )
    {
        p1 = (real *)in_data;
        p2 = (real *)out_data;
        
        /* n complex = 2*n real elements */
        for(i=0;i<2*n;i++)
        {
            p2[i] = p1[i];
        }
    }
  
    /* Elements 0   .. 2*n-1 in work are used for ffac values,
     * Elements 2*n .. 4*n-1 are internal FFTPACK work space.
     */
    
    if(dir == GMX_FFT_FORWARD) 
    {
        fftpack_cfftf1(n,(real *)out_data,fft->work+2*n,fft->work,fft->ifac, -1);
    }
    else if(dir == GMX_FFT_BACKWARD)
    {
        fftpack_cfftf1(n,(real *)out_data,fft->work+2*n,fft->work,fft->ifac, 1); 
    }
    else
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    

    return 0;
}



int 
gmx_fft_1d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           i,n;
    real *  p1;
    real *  p2;

    n = fft->n;
    
    if(n==1)
    {
        p1 = (real *)in_data;
        p2 = (real *)out_data;        
        p2[0] = p1[0];
        if(dir == GMX_FFT_REAL_TO_COMPLEX)
            p2[1] = 0.0;
    }
    
    if(dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        /* FFTPACK only does in-place transforms, so emulate out-of-place
         * by copying data to the output array first. This works fine, since
         * the complex array must be larger than the real.
         */
        if( in_data != out_data )
        {
            p1 = (real *)in_data;
            p2 = (real *)out_data;
            
            for(i=0;i<2*(n/2+1);i++)
            {
                p2[i] = p1[i];
            }
        }

        /* Elements 0 ..   n-1 in work are used for ffac values,
         * Elements n .. 2*n-1 are internal FFTPACK work space.
         */
        fftpack_rfftf1(n,(real *)out_data,fft->work+n,fft->work,fft->ifac);

        /*
         * FFTPACK has a slightly more compact storage than we, time to 
         * convert it: ove most of the array one step up to make room for 
         * zero imaginary parts. 
         */
        p2 = (real *)out_data;
        for(i=n-1;i>0;i--)
        {
            p2[i+1] = p2[i];
        }
        /* imaginary zero freq. */
        p2[1] = 0; 
        
        /* Is n even? */
        if( (n & 0x1) == 0 )
        {
            p2[n+1] = 0;
        }
        
    }
    else if(dir == GMX_FFT_COMPLEX_TO_REAL)
    {
        /* FFTPACK only does in-place transforms, and we cannot just copy
         * input to output first here since our real array is smaller than
         * the complex one. However, since the FFTPACK complex storage format 
         * is more compact than ours (2 reals) it will fit, so compact it
         * and copy on-the-fly to the output array.
         */
        p1 = (real *) in_data;
        p2 = (real *)out_data;

        p2[0] = p1[0];
        for(i=1;i<n;i++)
        {
            p2[i] = p1[i+1];
        }
        fftpack_rfftb1(n,(real *)out_data,fft->work+n,fft->work,fft->ifac); 
    }  
    else
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    
    
    return 0;
}


int 
gmx_fft_2d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int                i,nx,ny;
    t_complex *    data;
    
    nx = fft->n;
    ny = fft->next->n;
    
    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     * For 2D there is likely enough data to benefit from memcpy().
     */
    if( in_data != out_data )
    {
        memcpy(out_data,in_data,sizeof(t_complex)*nx*ny);
    }

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    /* y transforms */
    for(i=0;i<nx;i++)
    {
        gmx_fft_1d(fft->next,dir,data+i*ny,data+i*ny);
    }
    
    /* Transpose in-place to get data in place for x transform now */
    gmx_fft_transpose_2d(data,data,nx,ny);
    
    /* x transforms */
    for(i=0;i<ny;i++)
    {
        gmx_fft_1d(fft,dir,data+i*nx,data+i*nx);
    }
    
    /* Transpose in-place to get data back in original order */
    gmx_fft_transpose_2d(data,data,ny,nx);
    
    return 0;
}



int 
gmx_fft_2d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int                i,j,nx,ny,nyc;
    t_complex *    data;
    real *       work;
    real *       p1;
    real *       p2;
    
    nx=fft->n;
    ny=fft->next->n;
    /* Number of complex elements in y direction */
    nyc=(ny/2+1);
    
    work = fft->work+4*nx;
        
    if(dir==GMX_FFT_REAL_TO_COMPLEX)
    {
        /* If we are doing an in-place transform the 2D array is already
         * properly padded by the user, and we are all set.
         *
         * For out-of-place there is no array padding, but FFTPACK only
         * does in-place FFTs internally, so we need to start by copying
         * data from the input to the padded (larger) output array.
         */
        if( in_data != out_data )
        {
            p1 = (real *)in_data;
            p2 = (real *)out_data;
            
            for(i=0;i<nx;i++)
            {
                for(j=0;j<ny;j++)
                {
                    p2[i*nyc*2+j] = p1[i*ny+j];
                }
            }
        }
        data = (t_complex *)out_data;

        /* y real-to-complex FFTs */
        for(i=0;i<nx;i++)
        {
            gmx_fft_1d_real(fft->next,GMX_FFT_REAL_TO_COMPLEX,data+i*nyc,data+i*nyc);
        }
       
        /* Transform to get X data in place */
        gmx_fft_transpose_2d(data,data,nx,nyc);
        
        /* Complex-to-complex X FFTs */
        for(i=0;i<nyc;i++)
        {
            gmx_fft_1d(fft,GMX_FFT_FORWARD,data+i*nx,data+i*nx);
        }
  
        /* Transpose back */
        gmx_fft_transpose_2d(data,data,nyc,nx);
        
    }
    else if(dir==GMX_FFT_COMPLEX_TO_REAL)
    {
        /* An in-place complex-to-real transform is straightforward,
         * since the output array must be large enough for the padding to fit.
         *
         * For out-of-place complex-to-real transforms we cannot just copy
         * data to the output array, since it is smaller than the input.
         * In this case there's nothing to do but employing temporary work data,
         * starting at work+4*nx and using nx*nyc*2 elements.
         */
        if(in_data != out_data)
        {
            memcpy(work,in_data,sizeof(t_complex)*nx*nyc);
            data = (t_complex *)work;
        }
        else
        {
            /* in-place */
            data = (t_complex *)out_data;
        }

        /* Transpose to get X arrays */
        gmx_fft_transpose_2d(data,data,nx,nyc);
        
        /* Do X iFFTs */
        for(i=0;i<nyc;i++)
        {
            gmx_fft_1d(fft,GMX_FFT_BACKWARD,data+i*nx,data+i*nx);
        }
        
        /* Transpose to get Y arrays */
        gmx_fft_transpose_2d(data,data,nyc,nx);
        
        /* Do Y iFFTs */
        for(i=0;i<nx;i++)
        {
            gmx_fft_1d_real(fft->next,GMX_FFT_COMPLEX_TO_REAL,data+i*nyc,data+i*nyc);
        }
        
        if( in_data != out_data )
        {
            /* Output (pointed to by data) is now in padded format.
             * Pack it into out_data if we were doing an out-of-place transform.
             */
            p1 = (real *)data;
            p2 = (real *)out_data;
                
            for(i=0;i<nx;i++)
            {
                for(j=0;j<ny;j++)
                {
                    p2[i*ny+j] = p1[i*nyc*2+j];
                }
            }
        }
    }
    else
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    
    
    return 0;
}



int 
gmx_fft_3d          (gmx_fft_t                  fft,
                     enum gmx_fft_direction     dir,
                     void *                     in_data,
                     void *                     out_data)
{
    int              i,nx,ny,nz,rc;
    t_complex *  data;
    t_complex *  work;
    nx=fft->n;
    ny=fft->next->n;
    nz=fft->next->next->n;

    /* First 4*nx positions are FFTPACK workspace, then ours starts */
    work = (t_complex *)(fft->work+4*nx);
        
    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     * For 3D there is likely enough data to benefit from memcpy().
     */
    if( in_data != out_data )
    {
        memcpy(out_data,in_data,sizeof(t_complex)*nx*ny*nz);
    }
    
    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    /* Perform z transforms */
    for(i=0;i<nx*ny;i++)
        gmx_fft_1d(fft->next->next,dir,data+i*nz,data+i*nz);

    /* For each X slice, transpose the y & z dimensions inside the slice */
    for(i=0;i<nx;i++)
    {
        gmx_fft_transpose_2d(data+i*ny*nz,data+i*ny*nz,ny,nz);
    }

    /* Array is now (nx,nz,ny) - perform y transforms */
    for(i=0;i<nx*nz;i++)
    {
        gmx_fft_1d(fft->next,dir,data+i*ny,data+i*ny);
    }

    /* Transpose back to (nx,ny,nz) */
    for(i=0;i<nx;i++)
    {
        gmx_fft_transpose_2d(data+i*ny*nz,data+i*ny*nz,nz,ny);
    }
    
    /* Transpose entire x & y slices to go from
     * (nx,ny,nz) to (ny,nx,nz).
     * Use work data elements 4*n .. 4*n+2*nz-1. 
     */
    rc=gmx_fft_transpose_2d_nelem(data,data,nx,ny,nz,work);
    if( rc != 0)
    {
        gmx_fatal(FARGS,"Cannot transpose X & Y/Z in gmx_fft_3d().");
        return rc;
    }
    
    /* Then go from (ny,nx,nz) to (ny,nz,nx) */
    for(i=0;i<ny;i++)
    {
        gmx_fft_transpose_2d(data+i*nx*nz,data+i*nx*nz,nx,nz);
    }

    /* Perform x transforms */
    for(i=0;i<ny*nz;i++)
    {
        gmx_fft_1d(fft,dir,data+i*nx,data+i*nx);
    }
    
    /* Transpose back from (ny,nz,nx) to (ny,nx,nz) */
    for(i=0;i<ny;i++)
    {
        gmx_fft_transpose_2d(data+i*nz*nx,data+i*nz*nx,nz,nx);
    }
    
    /* Transpose from (ny,nx,nz) to (nx,ny,nz) 
     * Use work data elements 4*n .. 4*n+2*nz-1. 
     */
    rc = gmx_fft_transpose_2d_nelem(data,data,ny,nx,nz,work);
    if( rc != 0)
    {
        gmx_fatal(FARGS,"Cannot transpose Y/Z & X in gmx_fft_3d().");
        return rc;
    }
    
    return 0;
}


int 
gmx_fft_3d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int              i,j,k;
    int              nx,ny,nz,nzc,rc;
    t_complex *  data;
    t_complex *  work_transp;
    t_complex *  work_c2r;
    real *     p1;
    real *     p2;
    
    nx=fft->n;
    ny=fft->next->n;
    nz=fft->next->next->n;
    nzc=(nz/2+1);
        
    
    /* First 4*nx positions are FFTPACK workspace, then ours starts.
     * We have 2*nx*ny*nzc elements for temp complex-to-real storage when
     * doing out-of-place transforms, and another 2*nzc for transpose data.
     */
    work_c2r    = (t_complex *)(fft->work+4*nx);
    work_transp = (t_complex *)(fft->work+4*nx+2*nx*ny*nzc);

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    if(dir==GMX_FFT_REAL_TO_COMPLEX) 
    {
        /* FFTPACK only does in-place transforms, so emulate out-of-place
         * by copying data to the output array first. This is guaranteed to
         * work for real-to-complex since complex data is larger than the real.
         * For 3D there is likely enough data to benefit from memcpy().
         */
        if( in_data != out_data )
        {
            p1 = (real *)in_data;
            p2 = (real *)out_data;
           
            for(i=0;i<nx;i++)
            {
                for(j=0;j<ny;j++)
                {
                    for(k=0;k<nz;k++)
                    {
                        p2[(i*ny+j)*2*nzc+k] = p1[(i*ny+j)*nz+k];
                    }
                }
            }
        }
        data = (t_complex *)out_data;

        /* Transform the Y/Z slices real-to-complex */
        for(i=0;i<nx;i++)
        {
            gmx_fft_2d_real(fft->next,dir,data+i*ny*nzc,data+i*ny*nzc);
        }

        /* Transpose x & y slices to go from
         * (nx,ny,nzc) to (ny,nx,nzc).
         */
        rc=gmx_fft_transpose_2d_nelem(data,data,nx,ny,nzc,work_transp);
        if( rc != 0)
        {
            gmx_fatal(FARGS,"Cannot transpose X & Y/Z gmx_fft_3d_real().");
            return rc;
        }
        
        /* Then transpose from (ny,nx,nzc) to (ny,nzc,nx) */
        for(i=0;i<ny;i++)
        {
            gmx_fft_transpose_2d(data+i*nx*nzc,data+i*nx*nzc,nx,nzc);
        }
        
        /* Perform x transforms */
        for(i=0;i<ny*nzc;i++)
        {
            gmx_fft_1d(fft,GMX_FFT_FORWARD,data+i*nx,data+i*nx);
        }
        
        /* Transpose from (ny,nzc,nx) back to (ny,nx,nzc) */
        for(i=0;i<ny;i++)
        {
            gmx_fft_transpose_2d(data+i*nzc*nx,data+i*nzc*nx,nzc,nx);
        }
        
        /* Transpose back from (ny,nx,nzc) to (nx,ny,nz) */
        rc=gmx_fft_transpose_2d_nelem(data,data,ny,nx,nzc,work_transp);
        if( rc != 0)
        {
            gmx_fatal(FARGS,"Cannot transpose Y/Z & X in gmx_fft_3d_real().");
            return rc;
        }
            
    }
    else if(dir==GMX_FFT_COMPLEX_TO_REAL)
    {
        /* An in-place complex-to-real transform is straightforward,
         * since the output array must be large enough for the padding to fit.
         *
         * For out-of-place complex-to-real transforms we cannot just copy
         * data to the output array, since it is smaller than the input.
         * In this case there's nothing to do but employing temporary work data.
         */
        if(in_data != out_data)
        {
            memcpy(work_c2r,in_data,sizeof(t_complex)*nx*ny*nzc);
            data = (t_complex *)work_c2r;
        }
        else
        {
            /* in-place */
            data = (t_complex *)out_data;
        }
        
        /* Transpose x & y slices to go from
        * (nx,ny,nz) to (ny,nx,nz).
        */
        gmx_fft_transpose_2d_nelem(data,data,nx,ny,nzc,work_transp);

        /* Then go from (ny,nx,nzc) to (ny,nzc,nx) */
        for(i=0;i<ny;i++)
        {
            gmx_fft_transpose_2d(data+i*nx*nzc,data+i*nx*nzc,nx,nzc);
        }
        

        /* Perform x transforms */
        for(i=0;i<ny*nzc;i++)
        {
            gmx_fft_1d(fft,GMX_FFT_BACKWARD,data+i*nx,data+i*nx);
        }

        /* Transpose back from (ny,nzc,nx) to (ny,nx,nzc) */
        for(i=0;i<ny;i++)
        {
            gmx_fft_transpose_2d(data+i*nzc*nx,data+i*nzc*nx,nzc,nx);
        }
        
        /* Transpose back from (ny,nx,nzc) to (nx,ny,nz) */
        gmx_fft_transpose_2d_nelem(data,data,ny,nx,nzc,work_transp);
        

        /* Do 2D complex-to-real */
        for(i=0;i<nx;i++)
        {
            gmx_fft_2d_real(fft->next,dir,data+i*ny*nzc,data+i*ny*nzc);    
        }
        
        if( in_data != out_data )
        {
            /* Output (pointed to by data) is now in padded format.
             * Pack it into out_data if we were doing an out-of-place transform.
             */
            p1 = (real *)data;
            p2 = (real *)out_data;
            
            for(i=0;i<nx;i++)
            {
                for(j=0;j<ny;j++)
                {
                    for(k=0;k<nz;k++)
                    {
                        p2[(i*ny+j)*nz+k] = p1[(i*ny+j)*nzc*2+k];
                    }
                }
            }
        }
        
    }
    else
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    
    
    return 0;
}




void
gmx_fft_destroy(gmx_fft_t      fft)
{
    if(fft != NULL) 
    {
        free(fft->work);
        if(fft->next != NULL)
            gmx_fft_destroy(fft->next);
        free(fft);
    }
}
#else
int
gmx_fft_fftpack_empty;
#endif /* GMX_FFT_FFTPACK */
