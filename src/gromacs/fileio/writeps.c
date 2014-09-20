/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#include "gmxpre.h"

#include "writeps.h"

#include <stdio.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

const char *fontnm[efontNR] = {
    "Times-Roman", "Times-Italic",     "Times-Bold",    "Times-BoldItalic",
    "Helvetica",  "Helvetica-Oblique", "Helvetica-Bold", "Helvetica-BoldOblique",
    "Courier",    "Courier-Oblique",  "Courier-Bold",  "Courier-BoldOblique"
};


/* Internal psdata structure (abstract datatype)
 * to maintain the current state of the ps engine.
 */
struct t_int_psdata  {
    FILE   *fp;
    int     maxrgb;
    int     nrgb;
    t_rgb  *rgb;
    real    gen_ybox;
    int     ostack;
};


t_psdata ps_open(const char *fn, real x1, real y1, real x2, real y2)
{
    t_psdata ps;

    snew(ps, 1);

    ps->fp = gmx_fio_fopen(fn, "w");
    fprintf(ps->fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
    fprintf(ps->fp, "%%%%Creator: GROMACS\n");
    fprintf(ps->fp, "%%%%Title: %s\n", fn);
    fprintf(ps->fp, "%%%%BoundingBox: %g %g %g %g\n", x1, y1, x2, y2);
    fprintf(ps->fp, "%%%%EndComments\n");
    fprintf(ps->fp, "/m {moveto} bind def\n");
    fprintf(ps->fp, "/l {lineto} bind def\n");
    fprintf(ps->fp, "/rm {rmoveto} bind def\n");
    fprintf(ps->fp, "/r  {rlineto} bind def\n");
    fprintf(ps->fp, "/f {fill} bind def\n");
    fprintf(ps->fp, "/s {stroke} bind def\n");

    ps->nrgb     = 0;
    ps->maxrgb   = 0;
    ps->rgb      = NULL;
    ps->gen_ybox = 0;
    ps->ostack   = 0;

    return ps;
}

void ps_linewidth(t_psdata ps, int lw)
{
    fprintf(ps->fp, "%d setlinewidth\n", lw);
}

static void ps_defcolor(t_psdata ps, real r, real g, real b, char *cname)
{
    fprintf(ps->fp, "/%s {%g %g %g setrgbcolor} bind def\n", cname, r, g, b);
}

static void ps_selcolor(t_psdata ps, char *cname)
{
    fprintf(ps->fp, "%s\n", cname);
}

static int search_col(t_psdata ps, real r, real g, real b)
{
    int  i;
    char buf[12];

    for (i = 0; (i < ps->nrgb); i++)
    {
        if ((ps->rgb[i].r == r) && (ps->rgb[i].g == g) && (ps->rgb[i].b == b))
        {
            return i;
        }
    }

    if (ps->nrgb >= ps->maxrgb)
    {
        ps->maxrgb += 100;
        srenew(ps->rgb, ps->maxrgb);
    }

    sprintf(buf, "C%d", ps->nrgb);
    ps_defcolor(ps, r, g, b, buf);
    fprintf(ps->fp, "/B%d {%s b} bind def\n", ps->nrgb, buf);
    ps->rgb[i].r = r;
    ps->rgb[i].g = g;
    ps->rgb[i].b = b;
    ps->nrgb++;

    return (ps->nrgb-1);
}

void ps_color(t_psdata ps, real r, real g, real b)
{
    char buf[12];
    int  i;

    i = search_col(ps, r, g, b);

    sprintf(buf, "C%d", i);
    ps_selcolor(ps, buf);
}

void ps_rgb(t_psdata ps, t_rgb *rgb)
{
    ps_color(ps, rgb->r, rgb->g, rgb->b);
}


void ps_init_rgb_nbox(t_psdata ps, real xbox, real ybox)
{
    ps->gen_ybox = ybox;
    fprintf(ps->fp, "/by {def currentpoint "
            "%g y r %g %g r %g y neg r %g %g r f y add moveto} bind def\n",
            0.0, xbox, 0.0, 0.0, -xbox, 0.0);
    /* macro bn is used in ps_rgb_nbox to draw rectangular boxes */
}

void ps_rgb_nbox(t_psdata ps, t_rgb *rgb, real n)
{
    int i;

    if (n > 2)
    {
        ps_rgb(ps, rgb);
        fprintf(ps->fp, "/y %g by\n", n*ps->gen_ybox);
        /* macro by is defined in ps_init_rgb_nbox */
    }
    else
    {
        for (i = 0; (i < n); i++)
        {
            ps_rgb_box(ps, rgb);
        }
    }

}

void ps_init_rgb_box(t_psdata ps, real xbox, real ybox)
{
    fprintf(ps->fp, "/b {currentpoint "
            "%g %g r %g %g r %g %g r %g %g r f %g add moveto} bind def\n",
            0.0, ybox, xbox, 0.0, 0.0, -ybox, -xbox, 0.0, ybox);
    /* macro b is used in search_col to define macro B */
}

void ps_rgb_box(t_psdata ps, t_rgb *rgb)
{
    fprintf(ps->fp, "B%d\n", search_col(ps, rgb->r, rgb->g, rgb->b));
    /* macro B is defined in search_col from macro b */
}

void ps_lineto(t_psdata ps, real x, real y)
{
    fprintf(ps->fp, "%g %g l\n", x, y);
}

void ps_linerel(t_psdata ps, real dx, real dy)
{
    fprintf(ps->fp, "%g %g r\n", dx, dy);
}

void ps_moveto(t_psdata ps, real x, real y)
{
    fprintf(ps->fp, "%g %g m\n", x, y);
}

void ps_moverel(t_psdata ps, real dx, real dy)
{
    fprintf(ps->fp, "%g %g rm\n", dx, dy);
}

void ps_line(t_psdata ps, real x1, real y1, real x2, real y2)
{
    ps_moveto(ps, x1, y1);
    ps_lineto(ps, x2, y2);
    fprintf(ps->fp, "s\n");
}

static void do_box(t_psdata ps, real x1, real y1, real x2, real y2)
{
    ps_moveto(ps, x1, y1);
    ps_linerel(ps, 0, (real)(y2-y1));
    ps_linerel(ps, (real)(x2-x1), 0);
    ps_linerel(ps, 0, (real)(y1-y2));
    ps_linerel(ps, (real)(x1-x2), 0);
}

void ps_box(t_psdata ps, real x1, real y1, real x2, real y2)
{
    do_box(ps, x1, y1, x2, y2);
    fprintf(ps->fp, "s\n");
}

void ps_fillbox(t_psdata ps, real x1, real y1, real x2, real y2)
{
    do_box(ps, x1, y1, x2, y2);
    fprintf(ps->fp, "f\n");
}

void ps_arc(t_psdata ps, real x1, real y1, real rad, real a0, real a1)
{
    fprintf(ps->fp, "%g %g %g %g %g arc s\n", x1, y1, rad, a0, a1);
}

void ps_fillarc(t_psdata ps, real x1, real y1, real rad, real a0, real a1)
{
    fprintf(ps->fp, "%g %g %g %g %g arc f\n", x1, y1, rad, a0, a1);
}

void ps_arcslice(t_psdata ps, real xc, real yc,
                 real rad1, real rad2, real a0, real a1)
{
    fprintf(ps->fp, "newpath %g %g %g %g %g arc %g %g %g %g %g arcn closepath s\n",
            xc, yc, rad1, a0, a1, xc, yc, rad2, a1, a0);
}

void ps_fillarcslice(t_psdata ps, real xc, real yc,
                     real rad1, real rad2, real a0, real a1)
{
    fprintf(ps->fp, "newpath %g %g %g %g %g arc %g %g %g %g %g arcn closepath f\n",
            xc, yc, rad1, a0, a1, xc, yc, rad2, a1, a0);
}

void ps_circle(t_psdata ps, real x1, real y1, real rad)
{
    ps_arc(ps, x1, y1, rad, 0, 360);
}

void ps_font(t_psdata ps, int font, real size)
{

    if ((font < 0) || (font > efontNR))
    {
        fprintf(stderr, "Invalid Font: %d, using %s\n", font, fontnm[0]);
        font = 0;
    }
    fprintf(ps->fp, "/%s findfont\n", fontnm[font]);
    fprintf(ps->fp, "%g scalefont setfont\n", size);
}

void ps_strfont(t_psdata ps, char *font, real size)
{
    fprintf(ps->fp, "/%s findfont\n", font);
    fprintf(ps->fp, "%g scalefont setfont\n", size);
}

void ps_text(t_psdata ps, real x1, real y1, const char *str)
{
    ps_moveto(ps, x1, y1);
    fprintf(ps->fp, "(%s) show\n", str);
}

void ps_flip(t_psdata ps, gmx_bool bPlus)
{
    if (bPlus)
    {
        fprintf(ps->fp, "612.5 0 translate 90 rotate\n");
    }
    else
    {
        fprintf(ps->fp, "-90 rotate -612.5 0 translate\n");
    }
}

void ps_rotate(t_psdata ps, real angle)
{
    fprintf(ps->fp, "%f rotate\n", angle);
}

void ps_ctext(t_psdata ps, real x1, real y1, const char *str, int expos)
{
    if (expos == eXLeft)
    {
        ps_text(ps, x1, y1, str);
        return;
    }
    ps_moveto(ps, x1, y1);
    fprintf(ps->fp, "(%s) stringwidth\n", str);
    switch (expos)
    {
        case eXLeft:
            fprintf(ps->fp, "exch 0 exch pop exch\n");
            break;
        case eXCenter:
            fprintf(ps->fp, "exch 2 div neg exch\n");
            break;
        case eXRight:
            fprintf(ps->fp, "exch neg exch\n");
            break;
        default:
            gmx_fatal(FARGS, "invalid position index (expos=%d)", expos);
    }
    fprintf(ps->fp, "rmoveto (%s) show\n", str);
}

void ps_translate(t_psdata ps, real x, real y)
{
    fprintf(ps->fp, "%g %g translate\n", x, y);
}

void ps_setorigin(t_psdata ps)
{
    fprintf(ps->fp, "currentpoint dup 3 -1 roll dup 4 1 roll exch translate\n");
    ps->ostack++;
}

void ps_unsetorigin(t_psdata ps)
{
    if (ps->ostack <= 0)
    {
        gmx_fatal(FARGS, "No origin on stack!\n");
    }
    fprintf(ps->fp, "neg exch neg exch translate\n");
    ps->ostack--;
}

void ps_close(t_psdata ps)
{
    fprintf(ps->fp, "%%showpage\n");
    fprintf(ps->fp, "%%%%EOF\n");
    gmx_fio_fclose(ps->fp);
    sfree(ps->rgb);
    sfree(ps);
}

void ps_comment(t_psdata ps, const char *s)
{
    fprintf(ps->fp, "%%%% %s\n", s);
}
