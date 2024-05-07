/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <numeric>
#include <string>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/fileio/writeps.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#define FUDGE 1.2
#define DDD 2

typedef struct
{
    real     major;
    real     minor;
    real     offset;
    gmx_bool first;
    int      lineatzero;
    real     majorticklen;
    real     minorticklen;
    char     label[STRLEN];
    real     fontsize;
    char     font[STRLEN];
    real     tickfontsize;
    char     tickfont[STRLEN];
} t_axisdef;

typedef struct
{
    int       bw;
    real      linewidth;
    real      xoffs, yoffs;
    gmx_bool  bTitle;
    gmx_bool  bTitleOnce;
    gmx_bool  bYonce;
    real      titfontsize;
    char      titfont[STRLEN];
    gmx_bool  legend;
    real      legfontsize;
    char      legfont[STRLEN];
    char      leglabel[STRLEN];
    char      leg2label[STRLEN];
    real      xboxsize;
    real      yboxsize;
    real      boxspacing;
    real      boxlinewidth;
    real      ticklinewidth;
    real      zerolinewidth;
    t_axisdef X, Y;
} t_psrec;

/* MUST correspond to char *legend[] in main() */
enum
{
    elSel,
    elBoth,
    elFirst,
    elSecond,
    elNone,
    elNR
};

/* MUST correspond to char *combine[] in main() */
enum
{
    ecSel,
    ecHalves,
    ecAdd,
    ecSub,
    ecMult,
    ecDiv,
    ecNR
};

static void get_params(const char* mpin, const char* mpout, t_psrec* psr)
{
    static const char* gmx_bools[BOOL_NR + 1] = { "no", "yes", nullptr };
    /* this must correspond to t_rgb *linecolors[] below */
    static const char*     colors[] = { "none", "black", "white", nullptr };
    std::vector<t_inpfile> inp;

    WarningHandler wi{ false, 0 };

    if (mpin != nullptr)
    {
        std::string        libmpin = gmx::findLibraryFile(mpin).u8string();
        gmx::TextInputFile stream(libmpin);
        inp = read_inpfile(&stream, libmpin.c_str(), &wi);
    }
    else
    {
        inp.clear();
    }

    psr->bw        = get_eenum(&inp, "black&white", gmx_bools);
    psr->linewidth = get_ereal(&inp, "linewidth", 1.0, &wi);
    setStringEntry(&inp, "titlefont", psr->titfont, "Helvetica");
    psr->titfontsize = get_ereal(&inp, "titlefontsize", 20.0, &wi);
    psr->legend      = (get_eenum(&inp, "legend", gmx_bools) != 0);
    setStringEntry(&inp, "legendfont", psr->legfont, psr->titfont);
    setStringEntry(&inp, "legendlabel", psr->leglabel, "");
    setStringEntry(&inp, "legend2label", psr->leg2label, psr->leglabel);
    psr->legfontsize    = get_ereal(&inp, "legendfontsize", 14.0, &wi);
    psr->xboxsize       = get_ereal(&inp, "xbox", 0.0, &wi);
    psr->yboxsize       = get_ereal(&inp, "ybox", 0.0, &wi);
    psr->boxspacing     = get_ereal(&inp, "matrixspacing", 20.0, &wi);
    psr->xoffs          = get_ereal(&inp, "xoffset", 0.0, &wi);
    psr->yoffs          = get_ereal(&inp, "yoffset", psr->xoffs, &wi);
    psr->boxlinewidth   = get_ereal(&inp, "boxlinewidth", psr->linewidth, &wi);
    psr->ticklinewidth  = get_ereal(&inp, "ticklinewidth", psr->linewidth, &wi);
    psr->zerolinewidth  = get_ereal(&inp, "zerolinewidth", psr->ticklinewidth, &wi);
    psr->X.lineatzero   = get_eenum(&inp, "x-lineat0value", colors);
    psr->X.major        = get_ereal(&inp, "x-major", -1, &wi);
    psr->X.minor        = get_ereal(&inp, "x-minor", -1, &wi);
    psr->X.offset       = get_ereal(&inp, "x-firstmajor", 0.0, &wi);
    psr->X.first        = (get_eenum(&inp, "x-majorat0", gmx_bools) != 0);
    psr->X.majorticklen = get_ereal(&inp, "x-majorticklen", 8.0, &wi);
    psr->X.minorticklen = get_ereal(&inp, "x-minorticklen", 4.0, &wi);
    setStringEntry(&inp, "x-label", psr->X.label, "");
    psr->X.fontsize = get_ereal(&inp, "x-fontsize", 16.0, &wi);
    setStringEntry(&inp, "x-font", psr->X.font, psr->titfont);
    psr->X.tickfontsize = get_ereal(&inp, "x-tickfontsize", 10.0, &wi);
    setStringEntry(&inp, "x-tickfont", psr->X.tickfont, psr->X.font);
    psr->Y.lineatzero   = get_eenum(&inp, "y-lineat0value", colors);
    psr->Y.major        = get_ereal(&inp, "y-major", psr->X.major, &wi);
    psr->Y.minor        = get_ereal(&inp, "y-minor", psr->X.minor, &wi);
    psr->Y.offset       = get_ereal(&inp, "y-firstmajor", psr->X.offset, &wi);
    psr->Y.first        = (get_eenum(&inp, "y-majorat0", gmx_bools) != 0);
    psr->Y.majorticklen = get_ereal(&inp, "y-majorticklen", psr->X.majorticklen, &wi);
    psr->Y.minorticklen = get_ereal(&inp, "y-minorticklen", psr->X.minorticklen, &wi);
    setStringEntry(&inp, "y-label", psr->Y.label, psr->X.label);
    psr->Y.fontsize = get_ereal(&inp, "y-fontsize", psr->X.fontsize, &wi);
    setStringEntry(&inp, "y-font", psr->Y.font, psr->X.font);
    psr->Y.tickfontsize = get_ereal(&inp, "y-tickfontsize", psr->X.tickfontsize, &wi);
    setStringEntry(&inp, "y-tickfont", psr->Y.tickfont, psr->Y.font);

    check_warning_error(wi, FARGS);

    if (mpout != nullptr)
    {
        gmx::TextOutputFile stream(mpout);
        write_inpfile(&stream, mpout, &inp, TRUE, WriteMdpHeader::yes, &wi);
        stream.close();
    }

    done_warning(wi, FARGS);
}

static const t_rgb black = { 0, 0, 0 };
static const t_rgb white = { 1, 1, 1 };
#define BLACK (&black)
/* this must correspond to *colors[] in get_params */
static const t_rgb* const linecolors[] = { nullptr, &black, &white, nullptr };

static void leg_discrete(t_psdata*                      ps,
                         real                           x0,
                         real                           y0,
                         const std::string&             label,
                         real                           fontsize,
                         char*                          font,
                         gmx::ArrayRef<const t_mapping> map)
{
    real yhh;
    real boxhh;

    boxhh = fontsize + DDD;
    /* LANDSCAPE */
    ps_rgb(ps, BLACK);
    ps_strfont(ps, font, fontsize);
    yhh = y0 + fontsize + 3 * DDD;
    if (!label.empty())
    {
        ps_ctext(ps, x0, yhh, label, eXLeft);
    }
    ps_moveto(ps, x0, y0);
    for (const auto& m : map)
    {
        ps_setorigin(ps);
        ps_rgb(ps, &m.rgb);
        ps_fillbox(ps, DDD, DDD, DDD + fontsize, boxhh - DDD);
        ps_rgb(ps, BLACK);
        ps_box(ps, DDD, DDD, DDD + fontsize, boxhh - DDD);
        ps_ctext(ps, boxhh + 2 * DDD, fontsize / 3, m.desc, eXLeft);
        ps_unsetorigin(ps);
        ps_moverel(ps, DDD, -fontsize / 3);
    }
}

static void leg_continuous(t_psdata*                      ps,
                           real                           x0,
                           real                           x,
                           real                           y0,
                           const std::string&             label,
                           real                           fontsize,
                           char*                          font,
                           gmx::ArrayRef<const t_mapping> map,
                           int                            mapoffset)
{
    real       xx0;
    real       yhh, boxxh, boxyh;
    gmx::Index mapIndex = gmx::ssize(map) - mapoffset;

    boxyh = fontsize;
    if (x < 8 * fontsize)
    {
        x = 8 * fontsize;
    }
    boxxh = x / mapIndex;
    if (boxxh > fontsize)
    {
        boxxh = fontsize;
    }

    GMX_RELEASE_ASSERT(!map.empty(), "NULL map array provided to leg_continuous()");

    /* LANDSCAPE */
    xx0 = x0 - (mapIndex * boxxh) / 2.0;

    for (gmx::Index i = 0; (i < mapIndex); i++)
    {
        ps_rgb(ps, &(map[i + mapoffset].rgb));
        ps_fillbox(ps, xx0 + i * boxxh, y0, xx0 + (i + 1) * boxxh, y0 + boxyh);
    }
    ps_strfont(ps, font, fontsize);
    ps_rgb(ps, BLACK);
    ps_box(ps, xx0, y0, xx0 + mapIndex * boxxh, y0 + boxyh);

    yhh = y0 + boxyh + 3 * DDD;
    ps_ctext(ps, xx0 + boxxh / 2, yhh, map[0].desc, eXCenter);
    if (!label.empty())
    {
        ps_ctext(ps, x0, yhh, label, eXCenter);
    }
    ps_ctext(ps, xx0 + (mapIndex * boxxh) - boxxh / 2, yhh, map[map.size() - 1].desc, eXCenter);
}

static void leg_bicontinuous(t_psdata*                      ps,
                             real                           x0,
                             real                           x,
                             real                           y0,
                             const std::string&             label1,
                             const std::string&             label2,
                             real                           fontsize,
                             char*                          font,
                             gmx::ArrayRef<const t_mapping> map1,
                             gmx::ArrayRef<const t_mapping> map2)
{
    real xx1, xx2, x1, x2;

    x1  = x / (map1.size() + map2.size()) * map1.size(); /* width of legend 1 */
    x2  = x / (map1.size() + map2.size()) * map2.size(); /* width of legend 2 */
    xx1 = x0 - (x2 / 2.0) - fontsize;                    /* center of legend 1 */
    xx2 = x0 + (x1 / 2.0) + fontsize;                    /* center of legend 2 */
    x1 -= fontsize / 2;                                  /* adjust width */
    x2 -= fontsize / 2;                                  /* adjust width */
    leg_continuous(ps, xx1, x1, y0, label1, fontsize, font, map1, 0);
    leg_continuous(ps, xx2, x2, y0, label2, fontsize, font, map2, 0);
}

static real box_height(const t_matrix& mat, t_psrec* psr)
{
    return mat.ny * psr->yboxsize;
}

static real box_dh(t_psrec* psr)
{
    return psr->boxspacing;
}

static real box_dh_top(gmx_bool bOnce, t_psrec* psr)
{
    real dh;

    if (psr->bTitle || (psr->bTitleOnce && bOnce))
    {
        dh = 2 * psr->titfontsize;
    }
    else
    {
        dh = 0;
    }

    return dh;
}

static gmx_bool box_do_all_x_maj_ticks(t_psrec* psr)
{
    return (psr->boxspacing > (1.5 * psr->X.majorticklen));
}

static gmx_bool box_do_all_x_min_ticks(t_psrec* psr)
{
    return (psr->boxspacing > (1.5 * psr->X.minorticklen));
}

static void draw_boxes(t_psdata* ps, real x0, real y0, real w, gmx::ArrayRef<t_matrix> mat, t_psrec* psr)
{
    char   buf[128];
    real   xxx;
    char **xtick, **ytick;
    real   xx, yy, dy, xx00, yy00, offset_x, offset_y;
    int    x, ntx, nty;
    size_t strlength;

    /* Only necessary when there will be no y-labels */
    strlength = 0;

    /* Draw the box */
    ps_rgb(ps, BLACK);
    ps_linewidth(ps, static_cast<int>(psr->boxlinewidth));
    yy00 = y0;
    for (auto m = mat.begin(); m != mat.end(); ++m)
    {
        dy = box_height(*m, psr);
        ps_box(ps, x0 - 1, yy00 - 1, x0 + w + 1, yy00 + dy + 1);
        yy00 += dy + box_dh(psr) + box_dh_top(m + 1 == mat.end(), psr);
    }

    /* Draw the ticks on the axes */
    ps_linewidth(ps, static_cast<int>(psr->ticklinewidth));
    xx00         = x0 - 1;
    yy00         = y0 - 1;
    auto halfway = mat.begin() + (mat.size() / 2);
    for (auto m = mat.begin(); m != mat.end(); ++m)
    {
        if (m->flags & MAT_SPATIAL_X)
        {
            ntx      = m->nx + 1;
            offset_x = 0.1;
        }
        else
        {
            ntx      = m->nx;
            offset_x = 0.6;
        }
        if (m->flags & MAT_SPATIAL_Y)
        {
            nty      = m->ny + 1;
            offset_y = 0.1;
        }
        else
        {
            nty      = m->ny;
            offset_y = 0.6;
        }
        snew(xtick, ntx);
        for (int j = 0; (j < ntx); j++)
        {
            sprintf(buf, "%g", m->axis_x[j]);
            xtick[j] = gmx_strdup(buf);
        }
        ps_strfont(ps, psr->X.tickfont, psr->X.tickfontsize);
        for (x = 0; (x < ntx); x++)
        {
            xx = xx00 + (x + offset_x) * psr->xboxsize;
            if ((bRmod(m->axis_x[x], psr->X.offset, psr->X.major) || (psr->X.first && (x == 0)))
                && (m == mat.begin() || box_do_all_x_maj_ticks(psr)))
            {
                /* Longer tick marks */
                ps_line(ps, xx, yy00, xx, yy00 - psr->X.majorticklen);
                /* Plot label on lowest graph only */
                if (m == mat.begin())
                {
                    ps_ctext(ps, xx, yy00 - DDD - psr->X.majorticklen - psr->X.tickfontsize * 0.8, xtick[x], eXCenter);
                }
            }
            else if (bRmod(m->axis_x[x], psr->X.offset, psr->X.minor)
                     && ((m == mat.begin()) || box_do_all_x_min_ticks(psr)))
            {
                /* Shorter tick marks */
                ps_line(ps, xx, yy00, xx, yy00 - psr->X.minorticklen);
            }
            else if (bRmod(m->axis_x[x], psr->X.offset, psr->X.major))
            {
                /* Even shorter marks, only each X.major */
                ps_line(ps, xx, yy00, xx, yy00 - (psr->boxspacing / 2));
            }
        }
        ps_strfont(ps, psr->Y.tickfont, psr->Y.tickfontsize);
        snew(ytick, nty);
        for (int j = 0; (j < nty); j++)
        {
            sprintf(buf, "%g", m->axis_y[j]);
            ytick[j] = gmx_strdup(buf);
        }

        for (int y = 0; (y < nty); y++)
        {
            yy = yy00 + (y + offset_y) * psr->yboxsize;
            if (bRmod(m->axis_y[y], psr->Y.offset, psr->Y.major) || (psr->Y.first && (y == 0)))
            {
                /* Major ticks */
                strlength = std::max(strlength, std::strlen(ytick[y]));
                ps_line(ps, xx00, yy, xx00 - psr->Y.majorticklen, yy);
                ps_ctext(ps, xx00 - psr->Y.majorticklen - DDD, yy - psr->Y.tickfontsize / 3.0, ytick[y], eXRight);
            }
            else if (bRmod(m->axis_y[y], psr->Y.offset, psr->Y.minor))
            {
                /* Minor ticks */
                ps_line(ps, xx00, yy, xx00 - psr->Y.minorticklen, yy);
            }
        }
        sfree(xtick);
        sfree(ytick);

        /* Label on Y-axis */
        if (!psr->bYonce || m == halfway)
        {
            std::string mylab;
            if (strlen(psr->Y.label) > 0)
            {
                mylab = psr->Y.label;
            }
            else
            {
                mylab = m->label_y;
            }
            if (!mylab.empty())
            {
                ps_strfont(ps, psr->Y.font, psr->Y.fontsize);
                ps_flip(ps, TRUE);
                xxx = x0 - psr->X.majorticklen - psr->X.tickfontsize * strlength - DDD;
                ps_ctext(ps, yy00 + box_height(*m, psr) / 2.0, 612.5 - xxx, mylab, eXCenter);
                ps_flip(ps, FALSE);
            }
        }

        yy00 += box_height(*m, psr) + box_dh(psr) + box_dh_top(m + 1 == mat.end(), psr);
    }
    /* Label on X-axis */
    std::string mylab;
    if (strlen(psr->X.label) > 0)
    {
        mylab = psr->X.label;
    }
    else
    {
        mylab = mat[0].label_x;
    }
    if (!mylab.empty())
    {
        ps_strfont(ps, psr->X.font, psr->X.fontsize);
        ps_ctext(ps,
                 x0 + w / 2,
                 y0 - DDD - psr->X.majorticklen - psr->X.tickfontsize * FUDGE - psr->X.fontsize,
                 mylab,
                 eXCenter);
    }
}

static void draw_zerolines(t_psdata* out, real x0, real y0, real w, gmx::ArrayRef<t_matrix> mat, t_psrec* psr)
{
    real xx, yy, dy, xx00, yy00;
    int  x, y;

    xx00 = x0 - 1.5;
    yy00 = y0 - 1.5;
    ps_linewidth(out, static_cast<int>(psr->zerolinewidth));
    for (auto m = mat.begin(); m != mat.end(); ++m)
    {
        dy = box_height(*m, psr);
        /* m->axis_x and _y were already set by draw_boxes */
        if (psr->X.lineatzero)
        {
            ps_rgb(out, linecolors[psr->X.lineatzero]);
            for (x = 0; (x < m->nx); x++)
            {
                xx = xx00 + (x + 0.7) * psr->xboxsize;
                /* draw lines whenever tick label almost zero (e.g. next trajectory) */
                if (x != 0 && x < m->nx - 1
                    && std::abs(m->axis_x[x]) < 0.1 * std::abs(m->axis_x[x + 1] - m->axis_x[x]))
                {
                    ps_line(out, xx, yy00, xx, yy00 + dy + 2);
                }
            }
        }
        if (psr->Y.lineatzero)
        {
            ps_rgb(out, linecolors[psr->Y.lineatzero]);
            for (y = 0; (y < m->ny); y++)
            {
                yy = yy00 + (y + 0.7) * psr->yboxsize;
                /* draw lines whenever tick label almost zero (e.g. next trajectory) */
                if (y != 0 && y < m->ny - 1
                    && std::abs(m->axis_y[y]) < 0.1 * std::abs(m->axis_y[y + 1] - m->axis_y[y]))
                {
                    ps_line(out, xx00, yy, xx00 + w + 2, yy);
                }
            }
        }
        yy00 += box_height(*m, psr) + box_dh(psr) + box_dh_top(m + 1 == mat.end(), psr);
    }
}

static void box_dim(gmx::ArrayRef<t_matrix> mat,
                    gmx::ArrayRef<t_matrix> mat2,
                    t_psrec*                psr,
                    int                     elegend,
                    gmx_bool                bFrame,
                    real*                   w,
                    real*                   h,
                    real*                   dw,
                    real*                   dh)
{
    int  maxytick;
    real ww, hh, dww, dhh;

    hh = dww = dhh = 0;
    maxytick       = 0;

    ww = 0;
    for (const auto& m : mat)
    {
        ww = std::max(ww, m.nx * psr->xboxsize);
        hh += box_height(m, psr);
        maxytick = std::max(maxytick, m.nx);
    }
    if (bFrame)
    {
        if (mat[0].label_y[0])
        {
            dww += 2.0 * (psr->Y.fontsize + DDD);
        }
        if (psr->Y.major > 0)
        {
            dww += psr->Y.majorticklen + DDD
                   + psr->Y.tickfontsize * (std::log(static_cast<real>(maxytick)) / std::log(10.0));
        }
        else if (psr->Y.minor > 0)
        {
            dww += psr->Y.minorticklen;
        }

        if (mat[0].label_x[0])
        {
            dhh += psr->X.fontsize + 2 * DDD;
        }
        if (/* fool emacs auto-indent */
            (elegend == elBoth && (mat[0].legend[0] || (!mat2.empty() && mat2[0].legend[0])))
            || (elegend == elFirst && mat[0].legend[0])
            || (elegend == elSecond && (!mat2.empty() && mat2[0].legend[0])))
        {
            dhh += 2 * (psr->legfontsize * FUDGE + 2 * DDD);
        }
        else
        {
            dhh += psr->legfontsize * FUDGE + 2 * DDD;
        }
        if (psr->X.major > 0)
        {
            dhh += psr->X.tickfontsize * FUDGE + 2 * DDD + psr->X.majorticklen;
        }
        else if (psr->X.minor > 0)
        {
            dhh += psr->X.minorticklen;
        }

        hh += (mat.size() - 1) * box_dh(psr);
        hh += box_dh_top(TRUE, psr);
        if (mat.size() > 1)
        {
            hh += (mat.size() - 1) * box_dh_top(FALSE, psr);
        }
    }
    *w  = ww;
    *h  = hh;
    *dw = dww;
    *dh = dhh;
}

static std::vector<t_mapping> add_maps(gmx::ArrayRef<t_mapping> map1, gmx::ArrayRef<t_mapping> map2)
{
    static char mapper[] =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<"
            ".>/?";
    std::vector<t_mapping> map(map1.size() + map2.size());

    size_t nsymbols = std::strlen(mapper);
    if (map.size() > nsymbols * nsymbols)
    {
        gmx_fatal(FARGS, "Not enough symbols to merge the two colormaps\n");
    }
    printf("Combining colormaps of %zu and %zu elements into one of %zu elements\n",
           map1.size(),
           map2.size(),
           map.size());
    gmx::Index k = 0;
    for (gmx::Index j = 0; j < gmx::ssize(map1) && k < gmx::ssize(map); ++j, ++k)
    {
        map[k].code.c1 = mapper[k % nsymbols];
        if (map.size() > nsymbols)
        {
            map[k].code.c2 = mapper[k / nsymbols];
        }
        map[k].rgb.r = map1[j].rgb.r;
        map[k].rgb.g = map1[j].rgb.g;
        map[k].rgb.b = map1[j].rgb.b;
        map[k].desc  = map1[j].desc;
    }
    for (gmx::Index j = 0; j < gmx::ssize(map2) && k < gmx::ssize(map); ++j, ++k)
    {
        map[k].code.c1 = mapper[k % nsymbols];
        if (map.size() > nsymbols)
        {
            map[k].code.c2 = mapper[k / nsymbols];
        }
        map[k].rgb.r = map2[j].rgb.r;
        map[k].rgb.g = map2[j].rgb.g;
        map[k].rgb.b = map2[j].rgb.b;
        map[k].desc  = map2[j].desc;
    }

    return map;
}

static bool operator==(const t_mapping& lhs, const t_mapping& rhs)
{
    return (lhs.rgb.r == rhs.rgb.r && lhs.rgb.g == rhs.rgb.g && lhs.rgb.b == rhs.rgb.b);
}

static void xpm_mat(const char*             outf,
                    gmx::ArrayRef<t_matrix> mat,
                    gmx::ArrayRef<t_matrix> mat2,
                    gmx_bool                bDiag,
                    gmx_bool                bFirstDiag)
{
    FILE* out;
    int   x, y, col;

    out = gmx_ffopen(outf, "w");

    GMX_RELEASE_ASSERT(mat.size() == mat2.size(),
                       "Combined matrix write requires matrices of the same size");
    for (gmx::Index i = 0; i != gmx::ssize(mat); ++i)
    {
        // Color maps that differ only in RGB value are considered different
        if (mat2.empty() || std::equal(mat[i].map.begin(), mat[i].map.end(), mat2[i].map.begin()))
        {
            write_xpm_m(out, mat[0]);
        }
        else
        {
            auto map = add_maps(mat[i].map, mat2[i].map);
            for (x = 0; (x < mat[i].nx); x++)
            {
                for (y = 0; (y < mat[i].nx); y++)
                {
                    if ((x < y) || ((x == y) && bFirstDiag)) /* upper left  -> map1 */
                    {
                        col = mat[i].matrix(x, y);
                    }
                    else /* lower right -> map2 */
                    {
                        col = mat[i].map.size() + mat[i].matrix(x, y);
                    }
                    if ((bDiag) || (x != y))
                    {
                        mat[i].matrix(x, y) = col;
                    }
                    else
                    {
                        mat[i].matrix(x, y) = 0;
                    }
                }
            }
            mat[i].map = map;
            if (mat[i].title != mat2[i].title)
            {
                mat[i].title += " / " + mat2[i].title;
            }
            if (mat[i].legend != mat2[i].legend)
            {
                mat[i].legend += " / " + mat2[i].legend;
            }
            write_xpm_m(out, mat[i]);
        }
    }
    gmx_ffclose(out);
}

static void tick_spacing(int n, real axis[], real offset, char axisnm, real* major, real* minor)
{
    real     space;
    gmx_bool bTryAgain;
    int      i, j, t, f = 0, ten;
#define NFACT 4
    real major_fact[NFACT] = { 5, 4, 2, 1 };
    real minor_fact[NFACT] = { 5, 4, 4, 5 };

    /* start with interval between 10 matrix points: */
    space = std::max(10 * axis[1] - axis[0], axis[std::min(10, n - 1)] - axis[0]);
    /* get power of 10 */
    ten       = static_cast<int>(std::ceil(std::log(space) / std::log(10.0)) - 1);
    bTryAgain = TRUE;
    for (t = ten + 2; t > ten - 3 && bTryAgain; t--)
    {
        for (f = 0; f < NFACT && bTryAgain; f++)
        {
            space = std::pow(10.0_real, static_cast<real>(t)) * major_fact[f];
            /* count how many ticks we would get: */
            i = 0;
            for (j = 0; j < n; j++)
            {
                if (bRmod(axis[j], offset, space))
                {
                    i++;
                }
            }
            /* do we have a reasonable number of ticks ? */
            bTryAgain = (i > std::min(10, n - 1)) || (i < 5);
        }
    }
    if (bTryAgain)
    {
        space = std::max(10 * axis[1] - axis[0], axis[std::min(10, n - 1)] - axis[0]);
        fprintf(stderr, "Auto tick spacing failed for %c-axis, guessing %g\n", axisnm, space);
    }
    *major = space;
    *minor = space / minor_fact[(f > 0) ? f - 1 : 0];
    fprintf(stderr, "Auto tick spacing for %c-axis: major %g, minor %g\n", axisnm, *major, *minor);
}

static void ps_mat(const char*             outf,
                   gmx::ArrayRef<t_matrix> mat,
                   gmx::ArrayRef<t_matrix> mat2,
                   gmx_bool                bFrame,
                   gmx_bool                bDiag,
                   gmx_bool                bFirstDiag,
                   gmx_bool                bTitle,
                   gmx_bool                bTitleOnce,
                   gmx_bool                bYonce,
                   int                     elegend,
                   real                    size,
                   real                    boxx,
                   real                    boxy,
                   const char*             m2p,
                   const char*             m2pout,
                   int                     mapoffset)
{
    t_psrec  psrec, *psr;
    int      W, H;
    int      x, y, col;
    real     x0, y0, xx;
    real     w, h, dw, dh;
    gmx_bool bMap1, bNextMap1, bDiscrete;

    try
    {
        get_params(m2p, m2pout, &psrec);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    psr = &psrec;

    if (psr->X.major <= 0)
    {
        tick_spacing((mat[0].flags & MAT_SPATIAL_X) ? mat[0].nx + 1 : mat[0].nx,
                     mat[0].axis_x.data(),
                     psr->X.offset,
                     'X',
                     &(psr->X.major),
                     &(psr->X.minor));
    }
    if (psr->X.minor <= 0)
    {
        psr->X.minor = psr->X.major / 2;
    }
    if (psr->Y.major <= 0)
    {
        tick_spacing((mat[0].flags & MAT_SPATIAL_Y) ? mat[0].ny + 1 : mat[0].ny,
                     mat[0].axis_y.data(),
                     psr->Y.offset,
                     'Y',
                     &(psr->Y.major),
                     &(psr->Y.minor));
    }
    if (psr->Y.minor <= 0)
    {
        psr->Y.minor = psr->Y.major / 2;
    }

    if (boxx > 0)
    {
        psr->xboxsize = boxx;
        psr->yboxsize = boxx;
    }
    if (boxy > 0)
    {
        psr->yboxsize = boxy;
    }

    if (psr->xboxsize == 0)
    {
        psr->xboxsize = size / mat[0].nx;
        printf("Set the x-size of the box to %.3f\n", psr->xboxsize);
    }
    if (psr->yboxsize == 0)
    {
        psr->yboxsize = size / mat[0].nx;
        printf("Set the y-size of the box to %.3f\n", psr->yboxsize);
    }

    gmx::ArrayRef<const t_mapping> map1;
    int                            legendIndex = 0;
    for (const auto& m : mat)
    {
        if (m.map.size() > map1.size())
        {
            if (map1.empty())
            {
                printf("Selected legend of matrix # %d for display\n", legendIndex);
            }
            map1 = m.map;
        }
        ++legendIndex;
    }
    gmx::ArrayRef<const t_mapping> map2;
    if (!mat2.empty())
    {
        for (const auto& m : mat2)
        {
            if (m.map.size() > map2.size())
            {
                if (map2.empty())
                {
                    printf("Selected legend of matrix # %d for second display\n", legendIndex);
                }
                map2 = m.map;
            }
            ++legendIndex;
        }
    }
    if (mat[0].legend.empty() && psr->legend)
    {
        mat[0].legend = psr->leglabel;
    }

    bTitle          = bTitle && !mat.back().title.empty();
    bTitleOnce      = bTitleOnce && !mat.back().title.empty();
    psr->bTitle     = bTitle;
    psr->bTitleOnce = bTitleOnce;
    psr->bYonce     = bYonce;

    /* Set up size of box for nice colors */
    box_dim(mat, mat2, psr, elegend, bFrame, &w, &h, &dw, &dh);

    /* Set up bounding box */
    W = static_cast<int>(w + dw);
    H = static_cast<int>(h + dh);

    /* Start box at */
    x0 = dw;
    y0 = dh;
    x  = static_cast<int>(W + psr->xoffs);
    y  = static_cast<int>(H + psr->yoffs);
    if (bFrame)
    {
        x += 5 * DDD;
        y += 4 * DDD;
    }
    t_psdata out = ps_open(outf, 0, 0, x, y);
    ps_linewidth(&out, static_cast<int>(psr->linewidth));
    ps_init_rgb_box(&out, psr->xboxsize, psr->yboxsize);
    ps_init_rgb_nbox(&out, psr->xboxsize, psr->yboxsize);
    ps_translate(&out, psr->xoffs, psr->yoffs);

    if (bFrame)
    {
        ps_comment(&out, "Here starts the BOX drawing");
        draw_boxes(&out, x0, y0, w, mat, psr);
    }

    for (gmx::Index i = 0; i != gmx::ssize(mat); ++i)
    {
        if (bTitle || (bTitleOnce && i == gmx::ssize(mat) - 1))
        {
            /* Print title, if any */
            ps_rgb(&out, BLACK);
            ps_strfont(&out, psr->titfont, psr->titfontsize);
            std::string buf;
            if (mat2.empty() || mat[i].title == mat2[i].title)
            {
                buf = mat[i].title;
            }
            else
            {
                buf = mat[i].title + " / " + mat2[i].title;
            }
            ps_ctext(&out, x0 + w / 2, y0 + box_height(mat[i], psr) + psr->titfontsize, buf, eXCenter);
        }
        ps_comment(&out, gmx::formatString("Here starts the filling of box #%zd", i).c_str());
        for (x = 0; (x < mat[i].nx); x++)
        {
            int nexty;
            int nextcol;

            xx = x0 + x * psr->xboxsize;
            ps_moveto(&out, xx, y0);
            y     = 0;
            bMap1 = (mat2.empty() || (x < y || (x == y && bFirstDiag)));
            if ((bDiag) || (x != y))
            {
                col = mat[i].matrix(x, y);
            }
            else
            {
                col = -1;
            }
            for (nexty = 1; (nexty <= mat[i].ny); nexty++)
            {
                bNextMap1 = (mat2.empty() || (x < nexty || (x == nexty && bFirstDiag)));
                /* TRUE:  upper left  -> map1 */
                /* FALSE: lower right -> map2 */
                if ((nexty == mat[i].ny) || (!bDiag && (x == nexty)))
                {
                    nextcol = -1;
                }
                else
                {
                    nextcol = mat[i].matrix(x, nexty);
                }
                if ((nexty == mat[i].ny) || (col != nextcol) || (bMap1 != bNextMap1))
                {
                    if (col >= 0)
                    {
                        if (bMap1)
                        {
                            ps_rgb_nbox(&out, &(mat[i].map[col].rgb), nexty - y);
                        }
                        else
                        {
                            assert(!mat2.empty());
                            ps_rgb_nbox(&out, &(mat2[i].map[col].rgb), nexty - y);
                        }
                    }
                    else
                    {
                        ps_moverel(&out, 0, psr->yboxsize);
                    }
                    y     = nexty;
                    bMap1 = bNextMap1;
                    col   = nextcol;
                }
            }
        }
        y0 += box_height(mat[i], psr) + box_dh(psr) + box_dh_top(i + 1 == gmx::ssize(mat), psr);
    }

    if (psr->X.lineatzero || psr->Y.lineatzero)
    {
        /* reset y0 for first box */
        y0 = dh;
        ps_comment(&out, "Here starts the zero lines drawing");
        draw_zerolines(&out, x0, y0, w, mat, psr);
    }

    if (elegend != elNone)
    {
        std::string                    legend;
        gmx::ArrayRef<const t_mapping> leg_map;
        ps_comment(&out, "Now it's legend time!");
        ps_linewidth(&out, static_cast<int>(psr->linewidth));
        if (mat2.empty() || elegend != elSecond)
        {
            bDiscrete = mat[0].bDiscrete;
            legend    = mat[0].legend;
            leg_map   = map1;
        }
        else
        {
            bDiscrete = mat2[0].bDiscrete;
            legend    = mat2[0].legend;
            leg_map   = map2;
        }
        if (bDiscrete)
        {
            leg_discrete(&out, psr->legfontsize, DDD, legend, psr->legfontsize, psr->legfont, leg_map);
        }
        else
        {
            if (elegend != elBoth)
            {
                leg_continuous(
                        &out, x0 + w / 2, w / 2, DDD, legend, psr->legfontsize, psr->legfont, leg_map, mapoffset);
            }
            else
            {
                assert(!mat2.empty());
                leg_bicontinuous(
                        &out, x0 + w / 2, w, DDD, mat[0].legend, mat2[0].legend, psr->legfontsize, psr->legfont, map1, map2);
            }
        }
        ps_comment(&out, "Done processing");
    }

    ps_close(&out);
}

static void make_axis_labels(gmx::ArrayRef<t_matrix> mat1)
{
    for (auto& m : mat1)
    {
        /* Make labels for x axis */
        if (m.axis_x.empty())
        {
            m.axis_x.resize(m.nx);
            std::iota(m.axis_x.begin(), m.axis_x.end(), 0);
        }
        /* Make labels for y axis */
        if (m.axis_y.empty())
        {
            m.axis_y.resize(m.ny);
            std::iota(m.axis_y.begin(), m.axis_y.end(), 0);
        }
    }
}

static void prune_mat(gmx::ArrayRef<t_matrix> mat, gmx::ArrayRef<t_matrix> mat2, int skip)
{
    GMX_RELEASE_ASSERT(mat.size() == mat2.size() || mat2.empty(),
                       "Matrix pruning requires matrices of the same size");
    for (gmx::Index i = 0; i != gmx::ssize(mat); ++i)
    {
        fprintf(stderr,
                "converting %dx%d matrix to %dx%d\n",
                mat[i].nx,
                mat[i].ny,
                (mat[i].nx + skip - 1) / skip,
                (mat[i].ny + skip - 1) / skip);
        /* walk through matrix */
        int xs = 0;
        for (int x = 0; (x < mat[i].nx); x++)
        {
            if (x % skip == 0)
            {
                mat[i].axis_x[xs] = mat[i].axis_x[x];
                if (!mat2.empty())
                {
                    mat2[i].axis_x[xs] = mat2[i].axis_x[x];
                }
                int ys = 0;
                for (int y = 0; (y < mat[i].ny); y++)
                {
                    if (x == 0)
                    {
                        mat[i].axis_y[ys] = mat[i].axis_y[y];
                        if (!mat2.empty())
                        {
                            mat2[i].axis_y[ys] = mat2[i].axis_y[y];
                        }
                    }
                    if (y % skip == 0)
                    {
                        mat[i].matrix(xs, ys) = mat[i].matrix(x, y);
                        if (!mat2.empty())
                        {
                            mat2[i].matrix(xs, ys) = mat2[i].matrix(x, y);
                        }
                        ys++;
                    }
                }
                xs++;
            }
        }
        /* adjust parameters */
        mat[i].nx = (mat[i].nx + skip - 1) / skip;
        mat[i].ny = (mat[i].ny + skip - 1) / skip;
        if (!mat2.empty())
        {
            mat2[i].nx = (mat2[i].nx + skip - 1) / skip;
            mat2[i].ny = (mat2[i].ny + skip - 1) / skip;
        }
    }
}

static void zero_lines(gmx::ArrayRef<t_matrix> mat, gmx::ArrayRef<t_matrix> mat2)
{
    GMX_RELEASE_ASSERT(mat.size() == mat2.size(), "zero_lines requires matrices of the same size");
    for (gmx::Index i = 0; i != gmx::ssize(mat); ++i)
    {
        for (int m = 0; m < (!mat2.empty() ? 2 : 1); m++)
        {
            gmx::ArrayRef<t_matrix> mats;
            if (m == 0)
            {
                mats = mat;
            }
            else
            {
                mats = mat2;
            }
            for (int x = 0; x < mats[i].nx - 1; x++)
            {
                if (std::abs(mats[i].axis_x[x + 1]) < 1e-5)
                {
                    for (int y = 0; y < mats[i].ny; y++)
                    {
                        mats[i].matrix(x, y) = 0;
                    }
                }
            }
            for (int y = 0; y < mats[i].ny - 1; y++)
            {
                if (std::abs(mats[i].axis_y[y + 1]) < 1e-5)
                {
                    for (int x = 0; x < mats[i].nx; x++)
                    {
                        mats[i].matrix(x, y) = 0;
                    }
                }
            }
        }
    }
}

static void write_combined_matrix(int                     ecombine,
                                  const char*             fn,
                                  gmx::ArrayRef<t_matrix> mat1,
                                  gmx::ArrayRef<t_matrix> mat2,
                                  const real*             cmin,
                                  const real*             cmax)
{
    FILE*  out;
    real **rmat1, **rmat2;
    real   rhi, rlo;

    out = gmx_ffopen(fn, "w");
    GMX_RELEASE_ASSERT(mat1.size() == mat2.size(),
                       "Combined matrix write requires matrices of the same size");
    for (gmx::Index k = 0; k != gmx::ssize(mat1); k++)
    {
        if (mat2[k].nx != mat1[k].nx || mat2[k].ny != mat1[k].ny)
        {
            gmx_fatal(FARGS,
                      "Size of frame %zd in 1st (%dx%d) and 2nd matrix (%dx%d) do"
                      " not match.\n",
                      k,
                      mat1[k].nx,
                      mat1[k].ny,
                      mat2[k].nx,
                      mat2[k].ny);
        }
        printf("Combining two %dx%d matrices\n", mat1[k].nx, mat1[k].ny);
        rmat1 = matrix2real(&mat1[k], nullptr);
        rmat2 = matrix2real(&mat2[k], nullptr);
        if (nullptr == rmat1 || nullptr == rmat2)
        {
            gmx_fatal(FARGS,
                      "Could not extract real data from %s xpm matrices. Note that, e.g.,\n"
                      "gmx rms and gmx mdmat provide such data.\n",
                      (nullptr == rmat1 && nullptr == rmat2) ? "both" : "one of the");
        }
        rlo = 1e38;
        rhi = -1e38;
        for (int j = 0; j < mat1[k].ny; j++)
        {
            for (int i = 0; i < mat1[k].nx; i++)
            {
                switch (ecombine)
                {
                    case ecAdd: rmat1[i][j] += rmat2[i][j]; break;
                    case ecSub: rmat1[i][j] -= rmat2[i][j]; break;
                    case ecMult: rmat1[i][j] *= rmat2[i][j]; break;
                    case ecDiv: rmat1[i][j] /= rmat2[i][j]; break;
                    default: gmx_fatal(FARGS, "No such combination rule %d for matrices", ecombine);
                }
                rlo = std::min(rlo, rmat1[i][j]);
                rhi = std::max(rhi, rmat1[i][j]);
            }
        }
        if (cmin)
        {
            rlo = *cmin;
        }
        if (cmax)
        {
            rhi = *cmax;
        }
        int nlevels = static_cast<int>(std::max(mat1[k].map.size(), mat2[k].map.size()));
        if (rhi == rlo)
        {
            fprintf(stderr, "combination results in uniform matrix (%g), no output\n", rhi);
        }
        else
        {
            write_xpm(out,
                      mat1[k].flags,
                      mat1[k].title,
                      mat1[k].legend,
                      mat1[k].label_x,
                      mat1[k].label_y,
                      mat1[k].nx,
                      mat1[k].ny,
                      mat1[k].axis_x.data(),
                      mat1[k].axis_y.data(),
                      rmat1,
                      rlo,
                      rhi,
                      white,
                      black,
                      &nlevels);
        }
    }
    gmx_ffclose(out);
}

static void do_mat(gmx::ArrayRef<t_matrix> mat,
                   gmx::ArrayRef<t_matrix> mat2,
                   gmx_bool                bFrame,
                   gmx_bool                bZeroLine,
                   gmx_bool                bDiag,
                   gmx_bool                bFirstDiag,
                   gmx_bool                bTitle,
                   gmx_bool                bTitleOnce,
                   gmx_bool                bYonce,
                   int                     elegend,
                   real                    size,
                   real                    boxx,
                   real                    boxy,
                   const char*             epsfile,
                   const char*             xpmfile,
                   const char*             m2p,
                   const char*             m2pout,
                   int                     skip,
                   int                     mapoffset)
{
    GMX_RELEASE_ASSERT(mat.size() == mat2.size() || mat2.empty(),
                       "Combined matrix write requires matrices of the same size");
    if (!mat2.empty())
    {
        for (gmx::Index k = 0; k != gmx::ssize(mat); k++)
        {
            if ((mat2[k].nx != mat[k].nx) || (mat2[k].ny != mat[k].ny))
            {
                gmx_fatal(FARGS,
                          "WAKE UP!! Size of frame %zd in 2nd matrix file (%dx%d) does not match "
                          "size of 1st matrix (%dx%d) or the other way around.\n",
                          k,
                          mat2[k].nx,
                          mat2[k].ny,
                          mat[k].nx,
                          mat[k].ny);
            }
            for (int j = 0; (j < mat[k].ny); j++)
            {
                for (int i = bFirstDiag ? j + 1 : j; (i < mat[k].nx); i++)
                {
                    mat[k].matrix(i, j) = mat2[k].matrix(i, j);
                }
            }
        }
    }
    for (gmx::Index i = 0; i != gmx::ssize(mat); i++)
    {
        fprintf(stderr, "Matrix %zd is %d x %d\n", i, mat[i].nx, mat[i].ny);
    }

    make_axis_labels(mat);

    if (skip > 1)
    {
        prune_mat(mat, mat2, skip);
    }

    if (bZeroLine)
    {
        zero_lines(mat, mat);
    }

    if (epsfile != nullptr)
    {
        ps_mat(epsfile, mat, mat2, bFrame, bDiag, bFirstDiag, bTitle, bTitleOnce, bYonce, elegend, size, boxx, boxy, m2p, m2pout, mapoffset);
    }
    if (xpmfile != nullptr)
    {
        xpm_mat(xpmfile, mat, mat2, bDiag, bFirstDiag);
    }
}

static void gradient_map(const rvec grad, gmx::ArrayRef<t_mapping> map)
{
    int  i        = 0;
    real fraction = 1.0 / (map.size() - 1.0);
    for (auto& m : map)
    {
        real c  = i * fraction;
        m.rgb.r = 1 - c * (1 - grad[XX]);
        m.rgb.g = 1 - c * (1 - grad[YY]);
        m.rgb.b = 1 - c * (1 - grad[ZZ]);
        ++i;
    }
}

static void gradient_mat(rvec grad, gmx::ArrayRef<t_matrix> mat)
{
    for (auto& m : mat)
    {
        gradient_map(grad, m.map);
    }
}

static void rainbow_map(gmx_bool bBlue, gmx::ArrayRef<t_mapping> map)
{
    for (auto& m : map)
    {
        real c = (m.rgb.r + m.rgb.g + m.rgb.b) / 3;
        real r, g, b;
        if (c > 1)
        {
            c = 1;
        }
        if (bBlue)
        {
            c = 1 - c;
        }
        if (c <= 0.25) /* 0-0.25 */
        {
            r = 0;
            g = std::pow(4.0 * c, 2.0 / 3.0);
            b = 1;
        }
        else if (c <= 0.5) /* 0.25-0.5 */
        {
            r = 0;
            g = 1;
            b = std::pow(2.0 - 4.0 * c, 2.0 / 3.0);
        }
        else if (c <= 0.75) /* 0.5-0.75 */
        {
            r = std::pow(4.0 * c - 2.0, 2.0 / 3.0);
            g = 1;
            b = 0;
        }
        else /* 0.75-1 */
        {
            r = 1;
            g = std::pow(4.0 - 4.0 * c, 2.0 / 3.0);
            b = 0;
        }
        m.rgb.r = r;
        m.rgb.g = g;
        m.rgb.b = b;
    }
}

static void rainbow_mat(gmx_bool bBlue, gmx::ArrayRef<t_matrix> mat)
{
    for (auto& m : mat)
    {
        rainbow_map(bBlue, m.map);
    }
}

int gmx_xpm2ps(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] makes a beautiful color plot of an XPixelMap file.",
        "Labels and axis can be displayed, when they are supplied",
        "in the correct matrix format.",
        "Matrix data may be generated by programs such as [gmx-rms] or",
        "[gmx-mdmat].[PAR]",
        "Parameters are set in the [TT].m2p[tt] file optionally supplied with",
        "[TT]-di[tt]. Reasonable defaults are provided. Settings for the [IT]y[it]-axis",
        "default to those for the [IT]x[it]-axis. Font names have a defaulting hierarchy:",
        "titlefont -> legendfont; titlefont -> (xfont -> yfont -> ytickfont)",
        "-> xtickfont, e.g. setting titlefont sets all fonts, setting xfont",
        "sets yfont, ytickfont and xtickfont.[PAR]",
        "When no [TT].m2p[tt] file is supplied, many settings are taken from",
        "command line options. The most important option is [TT]-size[tt],",
        "which sets the size of the whole matrix in postscript units.",
        "This option can be overridden with the [TT]-bx[tt] and [TT]-by[tt]",
        "options (and the corresponding parameters in the [TT].m2p[tt] file),",
        "which set the size of a single matrix element.[PAR]",
        "With [TT]-f2[tt] a second matrix file can be supplied. Both matrix",
        "files will be read simultaneously and the upper left half of the",
        "first one ([TT]-f[tt]) is plotted together with the lower right",
        "half of the second one ([TT]-f2[tt]). The diagonal will contain",
        "values from the matrix file selected with [TT]-diag[tt].",
        "Plotting of the diagonal values can be suppressed altogether by",
        "setting [TT]-diag[tt] to [TT]none[tt].",
        "In this case, a new color map will be generated with",
        "a red gradient for negative numbers and a blue for positive.",
        "If the color coding and legend labels of both matrices are identical,",
        "only one legend will be displayed, else two separate legends are",
        "displayed.",
        "With [TT]-combine[tt], an alternative operation can be selected",
        "to combine the matrices. The output range is automatically set",
        "to the actual range of the combined matrix. This can be overridden",
        "with [TT]-cmin[tt] and [TT]-cmax[tt].[PAR]",
        "[TT]-title[tt] can be set to [TT]none[tt] to suppress the title, or to",
        "[TT]ylabel[tt] to show the title in the Y-label position (alongside",
        "the [IT]y[it]-axis).[PAR]",
        "With the [TT]-rainbow[tt] option, dull grayscale matrices can be turned",
        "into attractive color pictures.[PAR]",
        "Merged or rainbowed matrices can be written to an XPixelMap file with",
        "the [TT]-xpm[tt] option."
    };

    gmx_output_env_t* oenv;
    const char *      fn, *epsfile = nullptr, *xpmfile = nullptr;
    int               i, etitle, elegend, ediag, erainbow, ecombine;
    gmx_bool          bTitle, bTitleOnce, bDiag, bFirstDiag, bGrad;
    static gmx_bool   bFrame = TRUE, bZeroLine = FALSE, bYonce = FALSE;
    static real       size = 400, boxx = 0, boxy = 0, cmin = 0, cmax = 0;
    static rvec       grad = { 0, 0, 0 };
    enum
    {
        etSel,
        etTop,
        etOnce,
        etYlabel,
        etNone,
        etNR
    };
    const char* title[] = { nullptr, "top", "once", "ylabel", "none", nullptr };
    /* MUST correspond to enum elXxx as defined at top of file */
    const char* legend[] = { nullptr, "both", "first", "second", "none", nullptr };
    enum
    {
        edSel,
        edFirst,
        edSecond,
        edNone,
        edNR
    };
    const char* diag[] = { nullptr, "first", "second", "none", nullptr };
    enum
    {
        erSel,
        erNo,
        erBlue,
        erRed,
        erNR
    };
    const char* rainbow[] = { nullptr, "no", "blue", "red", nullptr };
    /* MUST correspond to enum ecXxx as defined at top of file */
    const char* combine[] = { nullptr, "halves", "add", "sub", "mult", "div", nullptr };
    static int  skip = 1, mapoffset = 0;
    t_pargs     pa[] = {
        { "-frame", FALSE, etBOOL, { &bFrame }, "Display frame, ticks, labels, title and legend" },
        { "-title", FALSE, etENUM, { title }, "Show title at" },
        { "-yonce", FALSE, etBOOL, { &bYonce }, "Show y-label only once" },
        { "-legend", FALSE, etENUM, { legend }, "Show legend" },
        { "-diag", FALSE, etENUM, { diag }, "Diagonal" },
        { "-size", FALSE, etREAL, { &size }, "Horizontal size of the matrix in ps units" },
        { "-bx",
          FALSE,
          etREAL,
          { &boxx },
          "Element x-size, overrides [TT]-size[tt] (also y-size when [TT]-by[tt] is not set)" },
        { "-by", FALSE, etREAL, { &boxy }, "Element y-size" },
        { "-rainbow", FALSE, etENUM, { rainbow }, "Rainbow colors, convert white to" },
        { "-gradient",
          FALSE,
          etRVEC,
          { grad },
          "Re-scale colormap to a smooth gradient from white {1,1,1} to {r,g,b}" },
        { "-skip", FALSE, etINT, { &skip }, "only write out every nr-th row and column" },
        { "-zeroline",
          FALSE,
          etBOOL,
          { &bZeroLine },
          "insert line in [REF].xpm[ref] matrix where axis label is zero" },
        { "-legoffset",
          FALSE,
          etINT,
          { &mapoffset },
          "Skip first N colors from [REF].xpm[ref] file for the legend" },
        { "-combine", FALSE, etENUM, { combine }, "Combine two matrices" },
        { "-cmin", FALSE, etREAL, { &cmin }, "Minimum for combination output" },
        { "-cmax", FALSE, etREAL, { &cmax }, "Maximum for combination output" }
    };
#define NPA asize(pa)
    t_filenm fnm[] = { { efXPM, "-f", nullptr, ffREAD },      { efXPM, "-f2", "root2", ffOPTRD },
                       { efM2P, "-di", nullptr, ffLIBOPTRD }, { efM2P, "-do", "out", ffOPTWR },
                       { efEPS, "-o", nullptr, ffOPTWR },     { efXPM, "-xpm", nullptr, ffOPTWR } };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW, NFILE, fnm, NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    etitle   = nenum(title);
    elegend  = nenum(legend);
    ediag    = nenum(diag);
    erainbow = nenum(rainbow);
    ecombine = nenum(combine);
    bGrad    = opt2parg_bSet("-gradient", NPA, pa);
    for (i = 0; i < DIM; i++)
    {
        if (grad[i] < 0 || grad[i] > 1)
        {
            gmx_fatal(FARGS, "RGB value %g out of range (0.0-1.0)", grad[i]);
        }
    }
    if (!bFrame)
    {
        etitle  = etNone;
        elegend = elNone;
    }

    epsfile = ftp2fn_null(efEPS, NFILE, fnm);
    xpmfile = opt2fn_null("-xpm", NFILE, fnm);
    if (epsfile == nullptr && xpmfile == nullptr)
    {
        if (ecombine != ecHalves)
        {
            xpmfile = opt2fn("-xpm", NFILE, fnm);
        }
        else
        {
            epsfile = ftp2fn(efEPS, NFILE, fnm);
        }
    }
    if (ecombine != ecHalves && epsfile)
    {
        fprintf(stderr,
                "WARNING: can only write result of arithmetic combination "
                "of two matrices to .xpm file\n"
                "         file %s will not be written\n",
                epsfile);
        epsfile = nullptr;
    }

    bDiag      = ediag != edNone;
    bFirstDiag = ediag != edSecond;

    fn = opt2fn("-f", NFILE, fnm);
    std::vector<t_matrix> mat, mat2;
    mat = read_xpm_matrix(fn);
    fprintf(stderr,
            "There %s %zu matri%s in %s\n",
            (mat.size() > 1) ? "are" : "is",
            mat.size(),
            (mat.size() > 1) ? "ces" : "x",
            fn);
    fn = opt2fn_null("-f2", NFILE, fnm);
    if (fn)
    {
        mat2 = read_xpm_matrix(fn);
        fprintf(stderr,
                "There %s %zu matri%s in %s\n",
                (mat2.size() > 1) ? "are" : "is",
                mat2.size(),
                (mat2.size() > 1) ? "ces" : "x",
                fn);
        if (mat.size() != mat2.size())
        {
            fprintf(stderr, "Different number of matrices, using the smallest number.\n");
            if (mat.size() > mat2.size())
            {
                mat.resize(mat2.size());
            }
            else
            {
                mat2.resize(mat.size());
            }
        }
    }
    else
    {
        if (ecombine != ecHalves)
        {
            fprintf(stderr,
                    "WARNING: arithmetic matrix combination selected (-combine), "
                    "but no second matrix (-f2) supplied\n"
                    "         no matrix combination will be performed\n");
        }
        ecombine = 0;
    }
    bTitle     = etitle == etTop;
    bTitleOnce = etitle == etOnce;
    if (etitle == etYlabel)
    {
        for (auto& m : mat)
        {
            m.label_y = m.title;
        }
        for (auto& m : mat2)
        {
            m.label_y = m.title;
        }
    }
    if (bGrad)
    {
        gradient_mat(grad, mat);
        if (!mat2.empty())
        {
            gradient_mat(grad, mat2);
        }
    }
    if (erainbow != erNo)
    {
        rainbow_mat(erainbow == erBlue, mat);
        if (!mat2.empty())
        {
            rainbow_mat(erainbow == erBlue, mat2);
        }
    }

    if (mat2.empty() && (elegend != elNone))
    {
        elegend = elFirst;
    }

    if (ecombine && ecombine != ecHalves)
    {
        write_combined_matrix(ecombine,
                              xpmfile,
                              mat,
                              mat2,
                              opt2parg_bSet("-cmin", NPA, pa) ? &cmin : nullptr,
                              opt2parg_bSet("-cmax", NPA, pa) ? &cmax : nullptr);
    }
    else
    {
        do_mat(mat,
               mat2,
               bFrame,
               bZeroLine,
               bDiag,
               bFirstDiag,
               bTitle,
               bTitleOnce,
               bYonce,
               elegend,
               size,
               boxx,
               boxy,
               epsfile,
               xpmfile,
               opt2fn_null("-di", NFILE, fnm),
               opt2fn_null("-do", NFILE, fnm),
               skip,
               mapoffset);
    }

    view_all(oenv, NFILE, fnm);

    return 0;
}
