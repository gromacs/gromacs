#include <cmath>

#include "../gmx_lapack.h"

#include "gromacs/utility/real.h"

#include "lapack_limits.h"

void F77_FUNC(sstebz, SSTEBZ)(const char* range,
                              const char* order,
                              int*        n,
                              float*      vl,
                              float*      vu,
                              int*        il,
                              int*        iu,
                              float*      abstol,
                              float*      d__,
                              float*      e,
                              int*        m,
                              int*        nsplit,
                              float*      w,
                              int*        iblock,
                              int*        isplit,
                              float*      work,
                              int*        iwork,
                              int*        info)
{
    int   i__1, i__2, i__3;
    float d__1, d__2, d__3, d__4, d__5;
    int   c__1 = 1;
    int   c__3 = 3;
    int   c__2 = 2;
    int   c__0 = 0;

    int         j, ib, jb, ie, je, nb;
    float       gl;
    int         im, in;
    float       gu;
    int         iw;
    float       wl, wu;
    int         nwl;
    float       ulp, wlu, wul;
    int         nwu;
    float       tmp1, tmp2;
    int         iend, ioff, iout, itmp1, jdisc;
    int         iinfo;
    float       atoli;
    int         iwoff;
    float       bnorm;
    int         itmax;
    float       wkill, rtoli, tnorm;
    int         ibegin;
    int         irange, idiscl;
    int         idumma[1];
    int         idiscu, iorder;
    int         ncnvrg;
    float       pivmin;
    int         toofew;
    const float safemn = GMX_FLOAT_MIN * (1.0 + GMX_FLOAT_EPS);

    --iwork;
    --work;
    --isplit;
    --iblock;
    --w;
    --e;
    --d__;

    *info = 0;

    if (*range == 'A' || *range == 'a')
    {
        irange = 1;
    }
    else if (*range == 'V' || *range == 'v')
    {
        irange = 2;
    }
    else if (*range == 'I' || *range == 'i')
    {
        irange = 3;
    }
    else
    {
        irange = 0;
    }

    if (*order == 'B' || *order == 'b')
    {
        iorder = 2;
    }
    else if (*order == 'E' || *order == 'e')
    {
        iorder = 1;
    }
    else
    {
        iorder = 0;
    }

    if (irange <= 0)
    {
        *info = -1;
    }
    else if (iorder <= 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (irange == 2)
    {
        if (*vl >= *vu)
        {
            *info = -5;
        }
    }
    else if (irange == 3 && (*il < 1 || *il > (*n)))
    {
        *info = -6;
    }
    else if (irange == 3 && (*iu < ((*n < *il) ? *n : *il) || *iu > *n))
    {
        *info = -7;
    }

    if (*info != 0)
    {
        return;
    }

    *info  = 0;
    ncnvrg = 0;
    toofew = 0;

    *m = 0;
    if (*n == 0)
    {
        return;
    }

    if (irange == 3 && *il == 1 && *iu == *n)
    {
        irange = 1;
    }

    ulp   = 2 * GMX_FLOAT_EPS;
    rtoli = ulp * 2.;
    nb    = DSTEBZ_BLOCKSIZE;
    if (nb <= 1)
    {
        nb = 0;
    }

    if (*n == 1)
    {
        *nsplit   = 1;
        isplit[1] = 1;
        if (irange == 2 && (*vl >= d__[1] || *vu < d__[1]))
        {
            *m = 0;
        }
        else
        {
            w[1]      = d__[1];
            iblock[1] = 1;
            *m        = 1;
        }
        return;
    }

    *nsplit  = 1;
    work[*n] = 0.;
    pivmin   = 1.;
    i__1     = *n;
    for (j = 2; j <= i__1; ++j)
    {
        d__1 = e[j - 1];
        tmp1 = d__1 * d__1;
        d__2 = ulp;
        if (std::abs(d__[j] * d__[j - 1]) * (d__2 * d__2) + safemn > tmp1)
        {
            isplit[*nsplit] = j - 1;
            ++(*nsplit);
            work[j - 1] = 0.;
        }
        else
        {
            work[j - 1] = tmp1;
            pivmin      = (pivmin > tmp1) ? pivmin : tmp1;
        }
    }
    isplit[*nsplit] = *n;
    pivmin *= safemn;

    if (irange == 3)
    {

        gu   = d__[1];
        gl   = d__[1];
        tmp1 = 0.;

        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j)
        {
            tmp2 = std::sqrt(work[j]);
            d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
            gu   = (d__1 > d__2) ? d__1 : d__2;
            d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
            gl   = (d__1 < d__2) ? d__1 : d__2;
            tmp1 = tmp2;
        }

        d__1 = gu, d__2 = d__[*n] + tmp1;
        gu   = (d__1 > d__2) ? d__1 : d__2;
        d__1 = gl, d__2 = d__[*n] - tmp1;
        gl    = (d__1 < d__2) ? d__1 : d__2;
        d__1  = std::abs(gl);
        d__2  = std::abs(gu);
        tnorm = (d__1 > d__2) ? d__1 : d__2;
        gl    = gl - tnorm * 2. * ulp * *n - pivmin * 4.;
        gu    = gu + tnorm * 2. * ulp * *n + pivmin * 2.;

        itmax = (int)((std::log(tnorm + pivmin) - std::log(pivmin)) / std::log(2.)) + 2;
        if (*abstol <= 0.)
        {
            atoli = ulp * tnorm;
        }
        else
        {
            atoli = *abstol;
        }

        work[*n + 1] = gl;
        work[*n + 2] = gl;
        work[*n + 3] = gu;
        work[*n + 4] = gu;
        work[*n + 5] = gl;
        work[*n + 6] = gu;
        iwork[1]     = -1;
        iwork[2]     = -1;
        iwork[3]     = *n + 1;
        iwork[4]     = *n + 1;
        iwork[5]     = *il - 1;
        iwork[6]     = *iu;

        F77_FUNC(slaebz, SLAEBZ)
        (&c__3,
         &itmax,
         n,
         &c__2,
         &c__2,
         &nb,
         &atoli,
         &rtoli,
         &pivmin,
         &d__[1],
         &e[1],
         &work[1],
         &iwork[5],
         &work[*n + 1],
         &work[*n + 5],
         &iout,
         &iwork[1],
         &w[1],
         &iblock[1],
         &iinfo);

        if (iwork[6] == *iu)
        {
            wl  = work[*n + 1];
            wlu = work[*n + 3];
            nwl = iwork[1];
            wu  = work[*n + 4];
            wul = work[*n + 2];
            nwu = iwork[4];
        }
        else
        {
            wl  = work[*n + 2];
            wlu = work[*n + 4];
            nwl = iwork[2];
            wu  = work[*n + 3];
            wul = work[*n + 1];
            nwu = iwork[3];
        }

        if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n)
        {
            *info = 4;
            return;
        }
    }
    else
    {


        /* avoid warnings for high gcc optimization */
        wlu = wul = 1.0;

        d__3  = std::abs(d__[1]) + std::abs(e[1]);
        d__4  = std::abs(d__[*n]) + std::abs(e[*n - 1]);
        tnorm = (d__3 > d__4) ? d__3 : d__4;

        i__1 = *n - 1;
        for (j = 2; j <= i__1; ++j)
        {
            d__4  = tnorm;
            d__5  = std::abs(d__[j]) + std::abs(e[j - 1]) + std::abs(e[j]);
            tnorm = (d__4 > d__5) ? d__4 : d__5;
        }

        if (*abstol <= 0.)
        {
            atoli = ulp * tnorm;
        }
        else
        {
            atoli = *abstol;
        }

        if (irange == 2)
        {
            wl = *vl;
            wu = *vu;
        }
        else
        {
            wl = 0.;
            wu = 0.;
        }
    }

    *m    = 0;
    iend  = 0;
    *info = 0;
    nwl   = 0;
    nwu   = 0;

    i__1 = *nsplit;
    for (jb = 1; jb <= i__1; ++jb)
    {
        ioff   = iend;
        ibegin = ioff + 1;
        iend   = isplit[jb];
        in     = iend - ioff;

        if (in == 1)
        {

            if (irange == 1 || wl >= d__[ibegin] - pivmin)
            {
                ++nwl;
            }
            if (irange == 1 || wu >= d__[ibegin] - pivmin)
            {
                ++nwu;
            }
            if (irange == 1 || ((wl < d__[ibegin] - pivmin) && (wu >= d__[ibegin] - pivmin)))
            {
                ++(*m);
                w[*m]      = d__[ibegin];
                iblock[*m] = jb;
            }
        }
        else
        {

            gu   = d__[ibegin];
            gl   = d__[ibegin];
            tmp1 = 0.;

            i__2 = iend - 1;
            for (j = ibegin; j <= i__2; ++j)
            {
                tmp2 = std::abs(e[j]);
                d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
                gu   = (d__1 > d__2) ? d__1 : d__2;
                d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
                gl   = (d__1 < d__2) ? d__1 : d__2;
                tmp1 = tmp2;
            }

            d__1 = gu, d__2 = d__[iend] + tmp1;
            gu   = (d__1 > d__2) ? d__1 : d__2;
            d__1 = gl, d__2 = d__[iend] - tmp1;
            gl    = (d__1 < d__2) ? d__1 : d__2;
            d__1  = std::abs(gl);
            d__2  = std::abs(gu);
            bnorm = (d__1 > d__2) ? d__1 : d__2;
            gl    = gl - bnorm * 2. * ulp * in - pivmin * 2.;
            gu    = gu + bnorm * 2. * ulp * in + pivmin * 2.;

            if (*abstol <= 0.)
            {
                d__1  = std::abs(gl);
                d__2  = std::abs(gu);
                atoli = ulp * ((d__1 > d__2) ? d__1 : d__2);
            }
            else
            {
                atoli = *abstol;
            }

            if (irange > 1)
            {
                if (gu < wl)
                {
                    nwl += in;
                    nwu += in;
                }
                gl = (gl > wl) ? gl : wl;
                gu = (gu < wu) ? gu : wu;
                if (gl >= gu) {}
                continue;
            }

            work[*n + 1]      = gl;
            work[*n + in + 1] = gu;
            F77_FUNC(slaebz, SLAEBZ)
            (&c__1,
             &c__0,
             &in,
             &in,
             &c__1,
             &nb,
             &atoli,
             &rtoli,
             &pivmin,
             &d__[ibegin],
             &e[ibegin],
             &work[ibegin],
             idumma,
             &work[*n + 1],
             &work[*n + (in << 1) + 1],
             &im,
             &iwork[1],
             &w[*m + 1],
             &iblock[*m + 1],
             &iinfo);

            nwl += iwork[1];
            nwu += iwork[in + 1];
            iwoff = *m - iwork[1];

            itmax = (int)((std::log(gu - gl + pivmin) - std::log(pivmin)) / std::log(2.)) + 2;
            F77_FUNC(slaebz, SLAEBZ)
            (&c__2,
             &itmax,
             &in,
             &in,
             &c__1,
             &nb,
             &atoli,
             &rtoli,
             &pivmin,
             &d__[ibegin],
             &e[ibegin],
             &work[ibegin],
             idumma,
             &work[*n + 1],
             &work[*n + (in << 1) + 1],
             &iout,
             &iwork[1],
             &w[*m + 1],
             &iblock[*m + 1],
             &iinfo);

            i__2 = iout;
            for (j = 1; j <= i__2; ++j)
            {
                tmp1 = (work[j + *n] + work[j + in + *n]) * .5;

                if (j > iout - iinfo)
                {
                    ncnvrg = 1;
                    ib     = -jb;
                }
                else
                {
                    ib = jb;
                }
                i__3 = iwork[j + in] + iwoff;
                for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je)
                {
                    w[je]      = tmp1;
                    iblock[je] = ib;
                }
            }

            *m += im;
        }
    }

    if (irange == 3)
    {
        im     = 0;
        idiscl = *il - 1 - nwl;
        idiscu = nwu - *iu;

        if (idiscl > 0 || idiscu > 0)
        {
            i__1 = *m;
            for (je = 1; je <= i__1; ++je)
            {
                if (w[je] <= wlu && idiscl > 0)
                {
                    --idiscl;
                }
                else if (w[je] >= wul && idiscu > 0)
                {
                    --idiscu;
                }
                else
                {
                    ++im;
                    w[im]      = w[je];
                    iblock[im] = iblock[je];
                }
            }
            *m = im;
        }
        if (idiscl > 0 || idiscu > 0)
        {

            if (idiscl > 0)
            {
                wkill = wu;
                i__1  = idiscl;
                for (jdisc = 1; jdisc <= i__1; ++jdisc)
                {
                    iw   = 0;
                    i__2 = *m;
                    for (je = 1; je <= i__2; ++je)
                    {
                        if (iblock[je] != 0 && (w[je] < wkill || iw == 0))
                        {
                            iw    = je;
                            wkill = w[je];
                        }
                    }
                    iblock[iw] = 0;
                }
            }
            if (idiscu > 0)
            {

                wkill = wl;
                i__1  = idiscu;
                for (jdisc = 1; jdisc <= i__1; ++jdisc)
                {
                    iw   = 0;
                    i__2 = *m;
                    for (je = 1; je <= i__2; ++je)
                    {
                        if (iblock[je] != 0 && (w[je] > wkill || iw == 0))
                        {
                            iw    = je;
                            wkill = w[je];
                        }
                    }
                    iblock[iw] = 0;
                }
            }
            im   = 0;
            i__1 = *m;
            for (je = 1; je <= i__1; ++je)
            {
                if (iblock[je] != 0)
                {
                    ++im;
                    w[im]      = w[je];
                    iblock[im] = iblock[je];
                }
            }
            *m = im;
        }
        if (idiscl < 0 || idiscu < 0)
        {
            toofew = 1;
        }
    }

    if (iorder == 1 && *nsplit > 1)
    {
        i__1 = *m - 1;
        for (je = 1; je <= i__1; ++je)
        {
            ie   = 0;
            tmp1 = w[je];
            i__2 = *m;
            for (j = je + 1; j <= i__2; ++j)
            {
                if (w[j] < tmp1)
                {
                    ie   = j;
                    tmp1 = w[j];
                }
            }

            if (ie != 0)
            {
                itmp1      = iblock[ie];
                w[ie]      = w[je];
                iblock[ie] = iblock[je];
                w[je]      = tmp1;
                iblock[je] = itmp1;
            }
        }
    }

    *info = 0;
    if (ncnvrg)
    {
        ++(*info);
    }
    if (toofew)
    {
        *info += 2;
    }
    return;
}
