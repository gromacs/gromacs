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
#ifndef GMX_FILEIO_WRITEPS_H
#define GMX_FILEIO_WRITEPS_H

#include <cstdio>

#include <string>
#include <vector>

#include "gromacs/fileio/rgb.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

typedef enum
{
    eXCenter,
    eXLeft,
    eXRight
} eXPos;

enum class Fonts : int
{
    Times,
    TimesItalic,
    TimesBold,
    TimesBoldItalic,
    Helvetica,
    HelveticaItalic,
    HelveticaBold,
    HelveticaBoldItalic,
    Courier,
    CourierItalic,
    CourierBold,
    CourierBoldItalic,
    Count
};


struct t_psdata
{
    FILE*              fp = nullptr;
    std::vector<t_rgb> rgb;
    real               gen_ybox = 0;
    int                ostack   = 0;
};


const char* enumValueToString(Fonts enumValue);

t_psdata ps_open(const char* fn, real x1, real y1, real x2, real y2);

void ps_linewidth(t_psdata* ps, int lw);
void ps_color(t_psdata* ps, real r, real g, real b);
void ps_rgb(t_psdata* ps, const t_rgb* rgb);

void ps_rgb_box(t_psdata* ps, t_rgb* rgb);
void ps_rgb_nbox(t_psdata* ps, t_rgb* rgb, real n);
void ps_init_rgb_box(t_psdata* ps, real xbox, real ybox);
void ps_init_rgb_nbox(t_psdata* ps, real xbox, real ybox);

void ps_lineto(t_psdata* ps, real x, real y);
void ps_linerel(t_psdata* ps, real dx, real dy);

void ps_moveto(t_psdata* ps, real x, real y);
void ps_moverel(t_psdata* ps, real dx, real dy);

void ps_line(t_psdata* ps, real x1, real y1, real x2, real y2);

void ps_box(t_psdata* ps, real x1, real y1, real x2, real y2);
void ps_fillbox(t_psdata* ps, real x1, real y1, real x2, real y2);

void ps_arc(t_psdata* ps, real x1, real y1, real rad, real a0, real a1);
void ps_fillarc(t_psdata* ps, real x1, real y1, real rad, real a0, real a1);
void ps_arcslice(t_psdata* ps, real xc, real yc, real rad1, real rad2, real a0, real a1);
void ps_fillarcslice(t_psdata* ps, real xc, real yc, real rad1, real rad2, real a0, real a1);

void ps_circle(t_psdata* ps, real x1, real y1, real rad);

void ps_font(t_psdata* ps, Fonts font, real size);
void ps_strfont(t_psdata* ps, char* font, real size);

void ps_text(t_psdata* ps, real x1, real y1, const std::string& str);
void ps_ctext(t_psdata* ps, real x1, real y1, const std::string& str, int expos);

void ps_close(t_psdata* ps);

void ps_flip(t_psdata* ps, gmx_bool bPlus);
/* Rotate over 90 (bPlus) or -90 (!bPlus) degrees */

void ps_rotate(t_psdata* ps, real angle);

void ps_translate(t_psdata* ps, real x, real y);

void ps_setorigin(t_psdata* ps);
void ps_unsetorigin(t_psdata* ps);

void ps_comment(t_psdata* ps, const char* s);

#endif
