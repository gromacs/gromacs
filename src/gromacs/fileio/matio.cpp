/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "matio.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"

static const char mapper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
#define NMAP (long int)strlen(mapper)

#define MAX_XPM_LINELENGTH 4096

real **mk_matrix(int nx, int ny, gmx_bool b1D)
{
    int    i;
    real **m;

    snew(m, nx);
    if (b1D)
    {
        snew(m[0], nx*ny);
    }

    for (i = 0; (i < nx); i++)
    {
        if (b1D)
        {
            m[i] = &(m[0][i*ny]);
        }
        else
        {
            snew(m[i], ny);
        }
    }

    return m;
}

void done_matrix(int nx, real ***m)
{
    int i;

    for (i = 0; (i < nx); i++)
    {
        sfree((*m)[i]);
    }
    sfree(*m);
    *m = nullptr;
}

gmx_bool matelmt_cmp(t_xpmelmt e1, t_xpmelmt e2)
{
    return (e1.c1 == e2.c1) && (e1.c2 == e2.c2);
}

t_matelmt searchcmap(int n, t_mapping map[], t_xpmelmt c)
{
    int i;

    for (i = 0; (i < n); i++)
    {
        if (matelmt_cmp(map[i].code, c))
        {
            return i;
        }
    }

    return -1;
}

int getcmap(FILE *in, const char *fn, t_mapping **map)
{
    int        i, n;
    char       line[STRLEN];
    char       code[STRLEN], desc[STRLEN];
    double     r, g, b;
    t_mapping *m;

    if (fgets2(line, STRLEN-1, in) == nullptr)
    {
        gmx_fatal(FARGS, "Not enough lines in colormap file %s"
                  "(just wanted to read number of entries)", fn);
    }
    sscanf(line, "%d", &n);
    snew(m, n);
    for (i = 0; (i < n); i++)
    {
        if (fgets2(line, STRLEN-1, in) == nullptr)
        {
            gmx_fatal(FARGS, "Not enough lines in colormap file %s"
                      "(should be %d, found only %d)", fn, n+1, i);
        }
        sscanf(line, "%s%s%lf%lf%lf", code, desc, &r, &g, &b);
        m[i].code.c1 = code[0];
        m[i].code.c2 = 0;
        m[i].desc    = gmx_strdup(desc);
        m[i].rgb.r   = r;
        m[i].rgb.g   = g;
        m[i].rgb.b   = b;
    }
    *map = m;

    return n;
}

int readcmap(const char *fn, t_mapping **map)
{
    FILE      *in;
    int        n;

    in = libopen(fn);
    n  = getcmap(in, fn, map);
    gmx_ffclose(in);

    return n;
}

void printcmap(FILE *out, int n, t_mapping map[])
{
    int i;

    fprintf(out, "%d\n", n);
    for (i = 0; (i < n); i++)
    {
        fprintf(out, "%c%c  %20s  %10g  %10g  %10g\n",
                map[i].code.c1 ? map[i].code.c1 : ' ',
                map[i].code.c2 ? map[i].code.c2 : ' ',
                map[i].desc, map[i].rgb.r, map[i].rgb.g, map[i].rgb.b);
    }
}

void writecmap(const char *fn, int n, t_mapping map[])
{
    FILE *out;

    out = gmx_fio_fopen(fn, "w");
    printcmap(out, n, map);
    gmx_fio_fclose(out);
}

static char *fgetline(char **line, int llmax, int *llalloc, FILE *in)
{
    char *fg;

    if (llmax > *llalloc)
    {
        srenew(*line, llmax+1);
        *llalloc = llmax;
    }
    fg = fgets(*line, llmax, in);
    trim(*line);

    return fg;
}

static void skipstr(char *line)
{
    int i, c;

    ltrim(line);
    c = 0;
    while ((line[c] != ' ') && (line[c] != '\0'))
    {
        c++;
    }
    i = c;
    while (line[c] != '\0')
    {
        line[c-i] = line[c];
        c++;
    }
    line[c-i] = '\0';
}

static char *line2string(char **line)
{
    int i;

    if (*line != nullptr)
    {
        while (((*line)[0] != '\"' ) && ( (*line)[0] != '\0' ))
        {
            (*line)++;
        }

        if ((*line)[0] != '\"')
        {
            return nullptr;
        }
        (*line)++;

        i = 0;
        while (( (*line)[i] != '\"' ) && ( (*line)[i] != '\0' ))
        {
            i++;
        }

        if ((*line)[i] != '\"')
        {
            *line = nullptr;
        }
        else
        {
            (*line)[i] = 0;
        }
    }

    return *line;
}

static void parsestring(char *line, const char *label, char *string)
{
    if (strstr(line, label))
    {
        if (strstr(line, label) < strchr(line, '\"'))
        {
            line2string(&line);
            strcpy(string, line);
        }
    }
}

static void read_xpm_entry(FILE *in, t_matrix *mm)
{
    t_mapping   *map;
    char        *line_buf = nullptr, *line = nullptr, *str, buf[256] = {0};
    int          i, m, col_len, nch = 0, n_axis_x, n_axis_y, llmax;
    int          llalloc = 0;
    unsigned int r, g, b;
    double       u;
    gmx_bool     bGetOnWithIt, bSetLine;
    t_xpmelmt    c;

    mm->flags      = 0;
    mm->title[0]   = 0;
    mm->legend[0]  = 0;
    mm->label_x[0] = 0;
    mm->label_y[0] = 0;
    mm->matrix     = nullptr;
    mm->axis_x     = nullptr;
    mm->axis_y     = nullptr;
    mm->bDiscrete  = FALSE;

    llmax = STRLEN;

    while ((nullptr != fgetline(&line_buf, llmax, &llalloc, in)) &&
           (std::strncmp(line_buf, "static", 6) != 0))
    {
        line = line_buf;
        parsestring(line, "title", (mm->title));
        parsestring(line, "legend", (mm->legend));
        parsestring(line, "x-label", (mm->label_x));
        parsestring(line, "y-label", (mm->label_y));
        parsestring(line, "type", buf);
    }

    if (!line || strncmp(line, "static", 6) != 0)
    {
        gmx_input("Invalid XPixMap");
    }

    if (buf[0] && (gmx_strcasecmp(buf, "Discrete") == 0))
    {
        mm->bDiscrete = TRUE;
    }

    if (debug)
    {
        fprintf(debug, "%s %s %s %s\n",
                mm->title, mm->legend, mm->label_x, mm->label_y);
    }

    /* Read sizes */
    bGetOnWithIt = FALSE;
    while (!bGetOnWithIt && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)))
    {
        line = line_buf;
        while (( line[0] != '\"' ) && ( line[0] != '\0' ))
        {
            line++;
        }

        if  (line[0] == '\"')
        {
            line2string(&line);
            sscanf(line, "%d %d %d %d", &(mm->nx), &(mm->ny), &(mm->nmap), &nch);
            if (nch > 2)
            {
                gmx_fatal(FARGS, "Sorry can only read xpm's with at most 2 caracters per pixel\n");
            }
            if (mm->nx <= 0 || mm->ny <= 0)
            {
                gmx_fatal(FARGS, "Dimensions of xpm-file have to be larger than 0\n");
            }
            llmax        = std::max(STRLEN, mm->nx+10);
            bGetOnWithIt = TRUE;
        }
    }
    if (debug)
    {
        fprintf(debug, "mm->nx %d mm->ny %d mm->nmap %d nch %d\n",
                mm->nx, mm->ny, mm->nmap, nch);
    }
    if (nch == 0)
    {
        gmx_fatal(FARGS, "Number of characters per pixel not found in xpm\n");
    }

    /* Read color map */
    snew(map, mm->nmap);
    m = 0;
    while ((m < mm->nmap) && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)))
    {
        line = std::strchr(line_buf, '\"');
        if  (line)
        {
            line++;
            /* Read xpm color map entry */
            map[m].code.c1 = line[0];
            if (nch == 1)
            {
                map[m].code.c2 = 0;
            }
            else
            {
                map[m].code.c2 = line[1];
            }
            line += nch;
            str   = std::strchr(line, '#');
            if (str)
            {
                str++;
                col_len = 0;
                while (std::isxdigit(str[col_len]))
                {
                    col_len++;
                }
                if (col_len == 6)
                {
                    sscanf(line, "%*s #%2x%2x%2x", &r, &g, &b);
                    map[m].rgb.r = r/255.0;
                    map[m].rgb.g = g/255.0;
                    map[m].rgb.b = b/255.0;
                }
                else if (col_len == 12)
                {
                    sscanf(line, "%*s #%4x%4x%4x", &r, &g, &b);
                    map[m].rgb.r = r/65535.0;
                    map[m].rgb.g = g/65535.0;
                    map[m].rgb.b = b/65535.0;
                }
                else
                {
                    gmx_file("Unsupported or invalid colormap in X PixMap");
                }
            }
            else
            {
                str = std::strchr(line, 'c');
                if (str)
                {
                    str += 2;
                }
                else
                {
                    gmx_file("Unsupported or invalid colormap in X PixMap");
                }
                fprintf(stderr, "Using white for color \"%s", str);
                map[m].rgb.r = 1;
                map[m].rgb.g = 1;
                map[m].rgb.b = 1;
            }
            line = std::strchr(line, '\"');
            line++;
            line2string(&line);
            map[m].desc = gmx_strdup(line);
            m++;
        }
    }
    if  (m != mm->nmap)
    {
        gmx_fatal(FARGS, "Number of read colors map entries (%d) does not match the number in the header (%d)", m, mm->nmap);
    }
    mm->map = map;

    /* Read axes, if there are any */
    n_axis_x     = 0;
    n_axis_y     = 0;
    bSetLine     = FALSE;
    do
    {
        if (bSetLine)
        {
            line = line_buf;
        }
        bSetLine = TRUE;
        if (strstr(line, "x-axis"))
        {
            line = std::strstr(line, "x-axis");
            skipstr(line);
            if (mm->axis_x == nullptr)
            {
                snew(mm->axis_x, mm->nx + 1);
            }
            while (sscanf(line, "%lf", &u) == 1)
            {
                if (n_axis_x > mm->nx)
                {
                    gmx_fatal(FARGS, "Too many x-axis labels in xpm (max %d)", mm->nx);
                }
                else if (n_axis_x == mm->nx)
                {
                    mm->flags |= MAT_SPATIAL_X;
                }
                mm->axis_x[n_axis_x] = u;
                n_axis_x++;
                skipstr(line);
            }
        }
        else if (std::strstr(line, "y-axis"))
        {
            line = std::strstr(line, "y-axis");
            skipstr(line);
            if (mm->axis_y == nullptr)
            {
                snew(mm->axis_y, mm->ny + 1);
            }
            while (sscanf(line, "%lf", &u) == 1)
            {
                if (n_axis_y > mm->ny)
                {
                    gmx_file("Too many y-axis labels in xpm");
                }
                else if (n_axis_y == mm->ny)
                {
                    mm->flags |= MAT_SPATIAL_Y;
                }
                mm->axis_y[n_axis_y] = u;
                n_axis_y++;
                skipstr(line);
            }
        }
    }
    while ((line[0] != '\"') && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)));

    /* Read matrix */
    snew(mm->matrix, mm->nx);
    for (i = 0; i < mm->nx; i++)
    {
        snew(mm->matrix[i], mm->ny);
    }
    m        = mm->ny-1;
    bSetLine = FALSE;
    do
    {
        if (bSetLine)
        {
            line = line_buf;
        }
        bSetLine = TRUE;
        if (m%(1+mm->ny/100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100*(mm->ny-m))/mm->ny);
        }
        while ((line[0] != '\"') && (line[0] != '\0'))
        {
            line++;
        }
        if (line[0] != '\"')
        {
            gmx_fatal(FARGS, "Not enough caracters in row %d of the matrix\n", m+1);
        }
        else
        {
            line++;
            for (i = 0; i < mm->nx; i++)
            {
                c.c1 = line[nch*i];
                if (nch == 1)
                {
                    c.c2 = 0;
                }
                else
                {
                    c.c2 = line[nch*i+1];
                }
                mm->matrix[i][m] = searchcmap(mm->nmap, mm->map, c);
            }
            m--;
        }
    }
    while ((m >= 0) && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)));
    if (m >= 0)
    {
        gmx_incons("Not enough rows in the matrix");
    }

    sfree(line_buf);
}

int read_xpm_matrix(const char *fnm, t_matrix **mat)
{
    FILE *in;
    char *line = nullptr;
    int   nmat;
    int   llalloc = 0;

    in = gmx_fio_fopen(fnm, "r");

    nmat = 0;
    while (nullptr != fgetline(&line, STRLEN, &llalloc, in))
    {
        if (std::strstr(line, "/* XPM */"))
        {
            srenew(*mat, nmat+1);
            read_xpm_entry(in, &(*mat)[nmat]);
            nmat++;
        }
    }
    gmx_fio_fclose(in);

    if (nmat == 0)
    {
        gmx_file("Invalid XPixMap");
    }

    sfree(line);

    return nmat;
}

real **matrix2real(t_matrix *in, real **out)
{
    t_mapping *map;
    double     tmp;
    real      *rmap;
    int        i, j, nmap;

    nmap = in->nmap;
    map  = in->map;
    snew(rmap, nmap);

    for (i = 0; i < nmap; i++)
    {
        if ((map[i].desc == nullptr) || (sscanf(map[i].desc, "%lf", &tmp) != 1))
        {
            fprintf(stderr, "Could not convert matrix to reals,\n"
                    "color map entry %d has a non-real description: \"%s\"\n",
                    i, map[i].desc);
            sfree(rmap);
            return nullptr;
        }
        rmap[i] = tmp;
    }

    if (out == nullptr)
    {
        snew(out, in->nx);
        for (i = 0; i < in->nx; i++)
        {
            snew(out[i], in->ny);
        }
    }
    for (i = 0; i < in->nx; i++)
    {
        for (j = 0; j < in->ny; j++)
        {
            out[i][j] = rmap[in->matrix[i][j]];
        }
    }

    sfree(rmap);

    fprintf(stderr, "Converted a %dx%d matrix with %d levels to reals\n",
            in->nx, in->ny, nmap);

    return out;
}

static void write_xpm_header(FILE *out,
                             const char *title, const char *legend,
                             const char *label_x, const char *label_y,
                             gmx_bool bDiscrete)
{
    fprintf(out,  "/* XPM */\n");
    try
    {
        gmx::BinaryInformationSettings settings;
        settings.generatedByHeader(true);
        settings.linePrefix("/* ");
        settings.lineSuffix(" */");
        gmx::printBinaryInformation(out, gmx::getProgramContext(),
                                    settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    fprintf(out,  "/* This file can be converted to EPS by the GROMACS program xpm2ps */\n");
    fprintf(out,  "/* title:   \"%s\" */\n", title);
    fprintf(out,  "/* legend:  \"%s\" */\n", legend);
    fprintf(out,  "/* x-label: \"%s\" */\n", label_x);
    fprintf(out,  "/* y-label: \"%s\" */\n", label_y);
    if (bDiscrete)
    {
        fprintf(out, "/* type:    \"Discrete\" */\n");
    }
    else
    {
        fprintf(out, "/* type:    \"Continuous\" */\n");
    }
}

static int calc_nmid(int nlevels, real lo, real mid, real hi)
{
    /* Take care that we have at least 1 entry in the mid to hi range
     */
    return std::min(std::max(0, static_cast<int>(((mid-lo)/(hi-lo))*(nlevels-1))),
                    nlevels-1);
}

static void write_xpm_map3(FILE *out, int n_x, int n_y, int *nlevels,
                           real lo, real mid, real hi,
                           t_rgb rlo, t_rgb rmid, t_rgb rhi)
{
    int    i, nmid;
    real   r, g, b, clev_lo, clev_hi;

    if (*nlevels > NMAP*NMAP)
    {
        fprintf(stderr, "Warning, too many levels (%d) in matrix, using %d only\n",
                *nlevels, (int)(NMAP*NMAP));
        *nlevels = NMAP*NMAP;
    }
    else if (*nlevels < 2)
    {
        fprintf(stderr, "Warning, too few levels (%d) in matrix, using 2 instead\n",
                *nlevels);
        *nlevels = 2;
    }
    if (!((mid >= lo) && (mid < hi)))
    {
        gmx_fatal(FARGS, "Lo: %f, Mid: %f, Hi: %f\n", lo, mid, hi);
    }

    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %d %d\",\n",
            n_x, n_y, *nlevels, (*nlevels <= NMAP) ? 1 : 2);

    nmid    = calc_nmid(*nlevels, lo, mid, hi);
    clev_lo = nmid;
    clev_hi = (*nlevels - 1 - nmid);
    for (i = 0; (i < nmid); i++)
    {
        r   = rlo.r+(i*(rmid.r-rlo.r)/clev_lo);
        g   = rlo.g+(i*(rmid.g-rlo.g)/clev_lo);
        b   = rlo.b+(i*(rmid.b-rlo.b)/clev_lo);
        fprintf(out, "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[i % NMAP],
                (*nlevels <= NMAP) ? ' ' : mapper[i/NMAP],
                static_cast<unsigned int>(round(255*r)),
                static_cast<unsigned int>(round(255*g)),
                static_cast<unsigned int>(round(255*b)),
                ((nmid - i)*lo + i*mid)/clev_lo);
    }
    for (i = 0; (i < (*nlevels-nmid)); i++)
    {
        r   = rmid.r+(i*(rhi.r-rmid.r)/clev_hi);
        g   = rmid.g+(i*(rhi.g-rmid.g)/clev_hi);
        b   = rmid.b+(i*(rhi.b-rmid.b)/clev_hi);
        fprintf(out, "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[(i+nmid) % NMAP],
                (*nlevels <= NMAP) ? ' ' : mapper[(i+nmid)/NMAP],
                static_cast<unsigned int>(round(255*r)),
                static_cast<unsigned int>(round(255*g)),
                static_cast<unsigned int>(round(255*b)),
                ((*nlevels - 1 - nmid - i)*mid + i*hi)/clev_hi);
    }
}

static void pr_simple_cmap(FILE *out, real lo, real hi, int nlevel, t_rgb rlo,
                           t_rgb rhi, int i0)
{
    int    i;
    real   r, g, b, fac;

    for (i = 0; (i < nlevel); i++)
    {
        fac = (i+1.0)/(nlevel);
        r   = rlo.r+fac*(rhi.r-rlo.r);
        g   = rlo.g+fac*(rhi.g-rlo.g);
        b   = rlo.b+fac*(rhi.b-rlo.b);
        fprintf(out, "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[(i+i0) % NMAP],
                (nlevel <= NMAP) ? ' ' : mapper[(i+i0)/NMAP],
                static_cast<unsigned int>(round(255*r)),
                static_cast<unsigned int>(round(255*g)),
                static_cast<unsigned int>(round(255*b)),
                lo+fac*(hi-lo));
    }
}

static void pr_discrete_cmap(FILE *out, int *nlevel, int i0)
{
    t_rgb  rgbd[16] = {
        { 1.0, 1.0, 1.0 }, /* white */
        { 1.0, 0.0, 0.0 }, /* red */
        { 1.0, 1.0, 0.0 }, /* yellow */
        { 0.0, 0.0, 1.0 }, /* blue */
        { 0.0, 1.0, 0.0 }, /* green */
        { 1.0, 0.0, 1.0 }, /* purple */
        { 1.0, 0.4, 0.0 }, /* orange */
        { 0.0, 1.0, 1.0 }, /* cyan */
        { 1.0, 0.4, 0.4 }, /* pink */
        { 1.0, 1.0, 0.0 }, /* yellow */
        { 0.4, 0.4, 1.0 }, /* lightblue */
        { 0.4, 1.0, 0.4 }, /* lightgreen */
        { 1.0, 0.4, 1.0 }, /* lightpurple */
        { 1.0, 0.7, 0.4 }, /* lightorange */
        { 0.4, 1.0, 1.0 }, /* lightcyan */
        { 0.0, 0.0, 0.0 }  /* black */
    };

    int    i, n;

    *nlevel = std::min(16, *nlevel);
    n       = *nlevel;
    for (i = 0; (i < n); i++)
    {
        fprintf(out, "\"%c%c c #%02X%02X%02X \" /* \"%3d\" */,\n",
                mapper[(i+i0) % NMAP],
                (n <= NMAP) ? ' ' : mapper[(i+i0)/NMAP],
                static_cast<unsigned int>(round(255*rgbd[i].r)),
                static_cast<unsigned int>(round(255*rgbd[i].g)),
                static_cast<unsigned int>(round(255*rgbd[i].b)),
                i);
    }
}



static void write_xpm_map_split(FILE *out, int n_x, int n_y,
                                int *nlevel_top, real lo_top, real hi_top,
                                t_rgb rlo_top, t_rgb rhi_top,
                                gmx_bool bDiscreteColor,
                                int *nlevel_bot, real lo_bot, real hi_bot,
                                t_rgb rlo_bot, t_rgb rhi_bot)
{
    int    ntot;

    ntot = *nlevel_top + *nlevel_bot;
    if (ntot > NMAP)
    {
        gmx_fatal(FARGS, "Warning, too many levels (%d) in matrix", ntot);
    }

    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %d %d\",\n", n_x, n_y, ntot, 1);

    if (bDiscreteColor)
    {
        pr_discrete_cmap(out, nlevel_bot, 0);
    }
    else
    {
        pr_simple_cmap(out, lo_bot, hi_bot, *nlevel_bot, rlo_bot, rhi_bot, 0);
    }

    pr_simple_cmap(out, lo_top, hi_top, *nlevel_top, rlo_top, rhi_top, *nlevel_bot);
}


static void write_xpm_map(FILE *out, int n_x, int n_y, int *nlevels,
                          real lo, real hi, t_rgb rlo, t_rgb rhi)
{
    int    i, nlo;
    real   invlevel, r, g, b;

    if (*nlevels > NMAP*NMAP)
    {
        fprintf(stderr, "Warning, too many levels (%d) in matrix, using %d only\n",
                *nlevels, (int)(NMAP*NMAP));
        *nlevels = NMAP*NMAP;
    }
    else if (*nlevels < 2)
    {
        fprintf(stderr, "Warning, too few levels (%d) in matrix, using 2 instead\n", *nlevels);
        *nlevels = 2;
    }

    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %d %d\",\n",
            n_x, n_y, *nlevels, (*nlevels <= NMAP) ? 1 : 2);

    invlevel = 1.0/(*nlevels-1);
    for (i = 0; (i < *nlevels); i++)
    {
        nlo = *nlevels-1-i;
        r   = (nlo*rlo.r+i*rhi.r)*invlevel;
        g   = (nlo*rlo.g+i*rhi.g)*invlevel;
        b   = (nlo*rlo.b+i*rhi.b)*invlevel;
        fprintf(out, "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[i % NMAP], (*nlevels <= NMAP) ? ' ' : mapper[i/NMAP],
                static_cast<unsigned int>(round(255*r)),
                static_cast<unsigned int>(round(255*g)),
                static_cast<unsigned int>(round(255*b)),
                (nlo*lo+i*hi)*invlevel);
    }
}

static void write_xpm_axis(FILE *out, const char *axis, gmx_bool bSpatial,
                           int n, real *label)
{
    int i;

    if (label)
    {
        for (i = 0; i < (bSpatial ? n+1 : n); i++)
        {
            if (i % 80 == 0)
            {
                if (i)
                {
                    fprintf(out, "*/\n");
                }
                fprintf(out, "/* %s-axis:  ", axis);
            }
            fprintf(out, "%g ", label[i]);
        }
        fprintf(out, "*/\n");
    }
}

static void write_xpm_data(FILE *out, int n_x, int n_y, real **mat,
                           real lo, real hi, int nlevels)
{
    int  i, j, c;
    real invlevel;

    invlevel = (nlevels-1)/(hi-lo);
    for (j = n_y-1; (j >= 0); j--)
    {
        if (j%(1+n_y/100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100*(n_y-j))/n_y);
        }
        fprintf(out, "\"");
        for (i = 0; (i < n_x); i++)
        {
            c = std::round((mat[i][j]-lo)*invlevel);
            if (c < 0)
            {
                c = 0;
            }
            if (c >= nlevels)
            {
                c = nlevels-1;
            }
            if (nlevels <= NMAP)
            {
                fprintf(out, "%c", mapper[c]);
            }
            else
            {
                fprintf(out, "%c%c", mapper[c % NMAP], mapper[c / NMAP]);
            }
        }
        if (j > 0)
        {
            fprintf(out, "\",\n");
        }
        else
        {
            fprintf(out, "\"\n");
        }
    }
}

static void write_xpm_data3(FILE *out, int n_x, int n_y, real **mat,
                            real lo, real mid, real hi, int nlevels)
{
    int  i, j, c = 0, nmid;
    real invlev_lo, invlev_hi;

    nmid      = calc_nmid(nlevels, lo, mid, hi);
    invlev_hi = (nlevels-1-nmid)/(hi-mid);
    invlev_lo = (nmid)/(mid-lo);

    for (j = n_y-1; (j >= 0); j--)
    {
        if (j%(1+n_y/100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100*(n_y-j))/n_y);
        }
        fprintf(out, "\"");
        for (i = 0; (i < n_x); i++)
        {
            if (mat[i][j] >= mid)
            {
                c = nmid+std::round((mat[i][j]-mid)*invlev_hi);
            }
            else if (mat[i][j] >= lo)
            {
                c = std::round((mat[i][j]-lo)*invlev_lo);
            }
            else
            {
                c = 0;
            }

            if (c < 0)
            {
                c = 0;
            }
            if (c >= nlevels)
            {
                c = nlevels-1;
            }
            if (nlevels <= NMAP)
            {
                fprintf(out, "%c", mapper[c]);
            }
            else
            {
                fprintf(out, "%c%c", mapper[c % NMAP], mapper[c / NMAP]);
            }
        }
        if (j > 0)
        {
            fprintf(out, "\",\n");
        }
        else
        {
            fprintf(out, "\"\n");
        }
    }
}

static void write_xpm_data_split(FILE *out, int n_x, int n_y, real **mat,
                                 real lo_top, real hi_top, int nlevel_top,
                                 real lo_bot, real hi_bot, int nlevel_bot)
{
    int  i, j, c;
    real invlev_top, invlev_bot;

    invlev_top = (nlevel_top-1)/(hi_top-lo_top);
    invlev_bot = (nlevel_bot-1)/(hi_bot-lo_bot);

    for (j = n_y-1; (j >= 0); j--)
    {
        if (j % (1+n_y/100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100*(n_y-j))/n_y);
        }
        fprintf(out, "\"");
        for (i = 0; (i < n_x); i++)
        {
            if (i < j)
            {
                c = nlevel_bot+round((mat[i][j]-lo_top)*invlev_top);
                if ((c < nlevel_bot) || (c >= nlevel_bot+nlevel_top))
                {
                    gmx_fatal(FARGS, "Range checking i = %d, j = %d, c = %d, bot = %d, top = %d matrix[i,j] = %f", i, j, c, nlevel_bot, nlevel_top, mat[i][j]);
                }
            }
            else if (i > j)
            {
                c = round((mat[i][j]-lo_bot)*invlev_bot);
                if ((c < 0) || (c >= nlevel_bot+nlevel_bot))
                {
                    gmx_fatal(FARGS, "Range checking i = %d, j = %d, c = %d, bot = %d, top = %d matrix[i,j] = %f", i, j, c, nlevel_bot, nlevel_top, mat[i][j]);
                }
            }
            else
            {
                c = nlevel_bot;
            }

            fprintf(out, "%c", mapper[c]);
        }
        if (j > 0)
        {
            fprintf(out, "\",\n");
        }
        else
        {
            fprintf(out, "\"\n");
        }
    }
}

void write_xpm_m(FILE *out, t_matrix m)
{
    /* Writes a t_matrix struct to .xpm file */

    int           i, j;
    gmx_bool      bOneChar;
    t_xpmelmt     c;

    bOneChar = (m.map[0].code.c2 == 0);
    write_xpm_header(out, m.title, m.legend, m.label_x, m.label_y,
                     m.bDiscrete);
    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %d %d\",\n", m.nx, m.ny, m.nmap, bOneChar ? 1 : 2);
    for (i = 0; (i < m.nmap); i++)
    {
        fprintf(out, "\"%c%c c #%02X%02X%02X \" /* \"%s\" */,\n",
                m.map[i].code.c1,
                bOneChar ? ' ' : m.map[i].code.c2,
                static_cast<unsigned int>(round(m.map[i].rgb.r*255)),
                static_cast<unsigned int>(round(m.map[i].rgb.g*255)),
                static_cast<unsigned int>(round(m.map[i].rgb.b*255)), m.map[i].desc);
    }
    write_xpm_axis(out, "x", m.flags & MAT_SPATIAL_X, m.nx, m.axis_x);
    write_xpm_axis(out, "y", m.flags & MAT_SPATIAL_Y, m.ny, m.axis_y);
    for (j = m.ny-1; (j >= 0); j--)
    {
        if (j%(1+m.ny/100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100*(m.ny-j))/m.ny);
        }
        fprintf(out, "\"");
        if (bOneChar)
        {
            for (i = 0; (i < m.nx); i++)
            {
                fprintf(out, "%c", m.map[m.matrix[i][j]].code.c1);
            }
        }
        else
        {
            for (i = 0; (i < m.nx); i++)
            {
                c = m.map[m.matrix[i][j]].code;
                fprintf(out, "%c%c", c.c1, c.c2);
            }
        }
        if (j > 0)
        {
            fprintf(out, "\",\n");
        }
        else
        {
            fprintf(out, "\"\n");
        }
    }
}

void write_xpm3(FILE *out, unsigned int flags,
                const std::string &title, const std::string &legend,
                const std::string &label_x, const std::string &label_y,
                int n_x, int n_y, real axis_x[], real axis_y[],
                real *mat[], real lo, real mid, real hi,
                t_rgb rlo, t_rgb rmid, t_rgb rhi, int *nlevels)
{
    /* See write_xpm.
     * Writes a colormap varying as rlo -> rmid -> rhi.
     */

    if (hi <= lo)
    {
        gmx_fatal(FARGS, "hi (%g) <= lo (%g)", hi, lo);
    }

    write_xpm_header(out, title.c_str(), legend.c_str(), label_x.c_str(), label_y.c_str(), FALSE);
    write_xpm_map3(out, n_x, n_y, nlevels, lo, mid, hi, rlo, rmid, rhi);
    write_xpm_axis(out, "x", flags & MAT_SPATIAL_X, n_x, axis_x);
    write_xpm_axis(out, "y", flags & MAT_SPATIAL_Y, n_y, axis_y);
    write_xpm_data3(out, n_x, n_y, mat, lo, mid, hi, *nlevels);
}

void write_xpm_split(FILE *out, unsigned int flags,
                     const std::string &title, const std::string &legend,
                     const std::string &label_x, const std::string &label_y,
                     int n_x, int n_y, real axis_x[], real axis_y[],
                     real *mat[],
                     real lo_top, real hi_top, int *nlevel_top,
                     t_rgb rlo_top, t_rgb rhi_top,
                     real lo_bot, real hi_bot, int *nlevel_bot,
                     gmx_bool bDiscreteColor,
                     t_rgb rlo_bot, t_rgb rhi_bot)
{
    /* See write_xpm.
     * Writes a colormap varying as rlo -> rmid -> rhi.
     */

    if (hi_top <= lo_top)
    {
        gmx_fatal(FARGS, "hi_top (%g) <= lo_top (%g)", hi_top, lo_top);
    }
    if (hi_bot <= lo_bot)
    {
        gmx_fatal(FARGS, "hi_bot (%g) <= lo_bot (%g)", hi_bot, lo_bot);
    }
    if (bDiscreteColor && (*nlevel_bot >= 16))
    {
        gmx_impl("Can not plot more than 16 discrete colors");
    }

    write_xpm_header(out, title.c_str(), legend.c_str(), label_x.c_str(), label_y.c_str(), FALSE);
    write_xpm_map_split(out, n_x, n_y, nlevel_top, lo_top, hi_top, rlo_top, rhi_top,
                        bDiscreteColor, nlevel_bot, lo_bot, hi_bot, rlo_bot, rhi_bot);
    write_xpm_axis(out, "x", flags & MAT_SPATIAL_X, n_x, axis_x);
    write_xpm_axis(out, "y", flags & MAT_SPATIAL_Y, n_y, axis_y);
    write_xpm_data_split(out, n_x, n_y, mat, lo_top, hi_top, *nlevel_top,
                         lo_bot, hi_bot, *nlevel_bot);
}

void write_xpm(FILE *out, unsigned int flags,
               const std::string &title, const std::string &legend,
               const std::string &label_x, const std::string &label_y,
               int n_x, int n_y, real axis_x[], real axis_y[],
               real *mat[], real lo, real hi,
               t_rgb rlo, t_rgb rhi, int *nlevels)
{
    /* out        xpm file
     * title      matrix title
     * legend     label for the continuous legend
     * label_x    label for the x-axis
     * label_y    label for the y-axis
     * n_x, n_y   size of the matrix
     * axis_x[]   the x-ticklabels
     * axis_y[]   the y-ticklables
     * *matrix[]  element x,y is matrix[x][y]
     * lo         output lower than lo is set to lo
     * hi         output higher than hi is set to hi
     * rlo        rgb value for level lo
     * rhi        rgb value for level hi
     * nlevels    number of color levels for the output
     */

    if (hi <= lo)
    {
        gmx_fatal(FARGS, "hi (%f) <= lo (%f)", hi, lo);
    }

    write_xpm_header(out, title.c_str(), legend.c_str(), label_x.c_str(), label_y.c_str(), FALSE);
    write_xpm_map(out, n_x, n_y, nlevels, lo, hi, rlo, rhi);
    write_xpm_axis(out, "x", flags & MAT_SPATIAL_X, n_x, axis_x);
    write_xpm_axis(out, "y", flags & MAT_SPATIAL_Y, n_y, axis_y);
    write_xpm_data(out, n_x, n_y, mat, lo, hi, *nlevels);
}
