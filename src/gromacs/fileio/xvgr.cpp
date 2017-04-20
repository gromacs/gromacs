/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "xvgr.h"

#include <cassert>
#include <cctype>
#include <cstring>

#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"

gmx_bool output_env_get_print_xvgr_codes(const gmx_output_env_t *oenv)
{
    int xvg_format;

    xvg_format = output_env_get_xvg_format(oenv);

    return (xvg_format == exvgXMGRACE || xvg_format == exvgXMGR);
}

static char *xvgrstr(const std::string &gmx, const gmx_output_env_t *oenv,
                     char *buf, int buflen)
{
    /* Supported greek letter names and corresponding xmgrace/xmgr symbols */
    const char *sym[]  = { "beta", "chi", "delta", "eta", "lambda", "mu", "omega", "phi", "psi", "rho", "theta", nullptr };
    const char  symc[] = { 'b',    'c',   'd',     'h',   'l',      'm',  'w',     'f',   'y',   'r',   'q',     '\0' };
    int         xvgf;
    gmx_bool    bXVGR;
    int         g, b, i;
    char        c;

    xvgf  = output_env_get_xvg_format(oenv);
    bXVGR = (xvgf == exvgXMGRACE || xvgf == exvgXMGR);

    g = 0;
    b = 0;
    while (gmx[g] != '\0')
    {
        /* Check with the largest string we have ("lambda"), add one for \0 */
        if (b + 6 + 1 >= buflen)
        {
            gmx_fatal(FARGS, "Output buffer length in xvgstr (%d) too small to process xvg input string '%s'", buflen, gmx.c_str());
        }
        if (gmx[g] == '\\')
        {
            g++;
            if (gmx[g] == 's')
            {
                /* Subscript */
                if (bXVGR)
                {
                    buf[b++] = '\\';
                    buf[b++] = 's';
                }
                else
                {
                    buf[b++] = '_';
                }
                g++;
            }
            else if (gmx[g] == 'S')
            {
                /* Superscript */
                if (bXVGR)
                {
                    buf[b++] = '\\';
                    buf[b++] = 'S';
                }
                else
                {
                    buf[b++] = '^';
                }
                g++;
            }
            else if (gmx[g] == 'N')
            {
                /* End sub/superscript */
                if (bXVGR)
                {
                    buf[b++] = '\\';
                    buf[b++] = 'N';
                }
                else
                {
                    if (gmx[g+1] != ' ')
                    {
                        buf[b++] = ' ';
                    }
                }
                g++;
            }
            else if (gmx[g] == '4')
            {
                /* Backward compatibility for xmgr normal font "\4" */
                switch (xvgf)
                {
                    case exvgXMGRACE:
                        sprintf(buf+b, "%s", "\\f{}");
                        break;
                    case exvgXMGR:
                        sprintf(buf+b, "%s", "\\4");
                        break;
                    default:
                        buf[b] = '\0';
                        break;
                }
                g++;
                b = strlen(buf);
            }
            else if (gmx[g] == '8')
            {
                /* Backward compatibility for xmgr symbol font "\8" */
                switch (xvgf)
                {
                    case exvgXMGRACE:
                        sprintf(buf+b, "%s", "\\x");
                        break;
                    case exvgXMGR:
                        sprintf(buf+b, "%s", "\\8");
                        break;
                    default:
                        buf[b] = '\0';
                        break;
                }
                g++;
                b = std::strlen(buf);
            }
            else
            {
                /* Check for special symbol */
                i = 0;
                while (sym[i] != nullptr &&
                       gmx_strncasecmp(sym[i], &gmx[g], std::strlen(sym[i])) != 0)
                {
                    i++;
                }
                if (sym[i] != nullptr)
                {
                    c = symc[i];
                    if (std::isupper(gmx[g]))
                    {
                        c = std::toupper(c);
                    }
                    switch (xvgf)
                    {
                        case exvgXMGRACE:
                            sprintf(buf+b, "%s%c%s", "\\x", c, "\\f{}");
                            break;
                        case exvgXMGR:
                            sprintf(buf+b, "%s%c%s", "\\8", c, "\\4");
                            break;
                        default:
                            std::strncat(buf+b, &gmx[g], std::strlen(sym[i]));
                            b += std::strlen(sym[i]);
                            if (gmx[g+std::strlen(sym[i])] != ' ')
                            {
                                buf[b++] = ' ';
                            }
                            buf[b] = '\0';
                            break;
                    }
                    g += std::strlen(sym[i]);
                    b  = std::strlen(buf);
                }
                else
                {
                    /* Unknown escape sequence, this should not happen.
                     * We do not generate a fatal error, since that might
                     * stop otherwise functioning code from working.
                     * Copy the backslash to the output and continue processing.
                     */
                    buf[b++] = '\\';
                }
            }
        }
        else
        {
            buf[b++] = gmx[g++];
        }
    }

    buf[b++] = '\0';

    return buf;
}

void xvgr_header(FILE *fp, const char *title, const std::string &xaxis,
                 const std::string &yaxis, int exvg_graph_type,
                 const gmx_output_env_t *oenv)
{
    char buf[STRLEN];

    if (output_env_get_print_xvgr_codes(oenv))
    {
        gmx_format_current_time(buf, STRLEN);
        fprintf(fp, "# This file was created %s", buf);
        try
        {
            gmx::BinaryInformationSettings settings;
            settings.generatedByHeader(true);
            settings.linePrefix("# ");
            gmx::printBinaryInformation(fp, output_env_get_program_context(oenv),
                                        settings);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        fprintf(fp, "# %s is part of G R O M A C S:\n#\n",
                output_env_get_program_display_name(oenv));
        fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
        fprintf(fp, "@    title \"%s\"\n", xvgrstr(title, oenv, buf, STRLEN));
        fprintf(fp, "@    xaxis  label \"%s\"\n",
                xvgrstr(xaxis, oenv, buf, STRLEN));
        fprintf(fp, "@    yaxis  label \"%s\"\n",
                xvgrstr(yaxis, oenv, buf, STRLEN));
        switch (exvg_graph_type)
        {
            case exvggtXNY:
                if (output_env_get_xvg_format(oenv) == exvgXMGR)
                {
                    fprintf(fp, "@TYPE nxy\n");
                }
                else
                {
                    fprintf(fp, "@TYPE xy\n");
                }
                break;
            case exvggtXYDY:
                fprintf(fp, "@TYPE xydy\n");
                break;
            case exvggtXYDYDY:
                fprintf(fp, "@TYPE xydydy\n");
                break;
        }
    }
}

FILE *xvgropen_type(const char *fn, const char *title, const std::string &xaxis,
                    const std::string &yaxis, int exvg_graph_type,
                    const gmx_output_env_t *oenv)
{
    FILE  *fp;

    fp = gmx_fio_fopen(fn, "w");

    xvgr_header(fp, title, xaxis, yaxis, exvg_graph_type, oenv);

    return fp;
}

FILE *xvgropen(const char *fn, const char *title, const std::string &xaxis,
               const std::string &yaxis, const gmx_output_env_t *oenv)
{
    return xvgropen_type(fn, title, xaxis, yaxis, exvggtXNY, oenv);
}

void
xvgrclose(FILE *fp)
{
    gmx_fio_fclose(fp);
}

void xvgr_subtitle(FILE *out, const char *subtitle, const gmx_output_env_t *oenv)
{
    char buf[STRLEN];

    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@ subtitle \"%s\"\n", xvgrstr(subtitle, oenv, buf, STRLEN));
    }
}

void xvgr_view(FILE *out, real xmin, real ymin, real xmax, real ymax,
               const gmx_output_env_t *oenv)
{
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@ view %g, %g, %g, %g\n", xmin, ymin, xmax, ymax);
    }
}

void xvgr_world(FILE *out, real xmin, real ymin, real xmax, real ymax,
                const gmx_output_env_t *oenv)
{
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@ world xmin %g\n"
                "@ world ymin %g\n"
                "@ world xmax %g\n"
                "@ world ymax %g\n", xmin, ymin, xmax, ymax);
    }
}

void xvgr_legend(FILE *out, int nsets, const char** setname,
                 const gmx_output_env_t *oenv)
{
    int  i;
    char buf[STRLEN];

    if (output_env_get_print_xvgr_codes(oenv))
    {
        xvgr_view(out, 0.15, 0.15, 0.75, 0.85, oenv);
        fprintf(out, "@ legend on\n");
        fprintf(out, "@ legend box on\n");
        fprintf(out, "@ legend loctype view\n");
        fprintf(out, "@ legend %g, %g\n", 0.78, 0.8);
        fprintf(out, "@ legend length %d\n", 2);
        for (i = 0; (i < nsets); i++)
        {
            if (setname[i])
            {
                if (output_env_get_xvg_format(oenv) == exvgXMGR)
                {
                    fprintf(out, "@ legend string %d \"%s\"\n",
                            i, xvgrstr(setname[i], oenv, buf, STRLEN));
                }
                else
                {
                    fprintf(out, "@ s%d legend \"%s\"\n",
                            i, xvgrstr(setname[i], oenv, buf, STRLEN));
                }
            }
        }
    }
}

void xvgr_new_dataset(FILE *out, int nr_first, int nsets,
                      const char **setname,
                      const gmx_output_env_t *oenv)
{
    int  i;
    char buf[STRLEN];

    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@\n");
        for (i = 0; (i < nsets); i++)
        {
            if (setname[i])
            {
                if (output_env_get_xvg_format(oenv) == exvgXMGR)
                {
                    fprintf(out, "@ legend string %d \"%s\"\n",
                            i+nr_first, xvgrstr(setname[i], oenv, buf, STRLEN));
                }
                else
                {
                    fprintf(out, "@ s%d legend \"%s\"\n",
                            i+nr_first, xvgrstr(setname[i], oenv, buf, STRLEN));
                }
            }
        }
    }
    else
    {
        fprintf(out, "\n");
    }
}

void xvgr_line_props(FILE *out, int NrSet, int LineStyle, int LineColor,
                     const gmx_output_env_t *oenv)
{
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@    with g0\n");
        fprintf(out, "@    s%d linestyle %d\n", NrSet, LineStyle);
        fprintf(out, "@    s%d color %d\n", NrSet, LineColor);
    }
}

static const char *LocTypeStr[] = { "view", "world" };
static const char *BoxFillStr[] = { "none", "color", "pattern" };

void xvgr_box(FILE *out,
              int LocType,
              real xmin, real ymin, real xmax, real ymax,
              int LineStyle, int LineWidth, int LineColor,
              int BoxFill, int BoxColor, int BoxPattern, const gmx_output_env_t *oenv)
{
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@with box\n");
        fprintf(out, "@    box on\n");
        fprintf(out, "@    box loctype %s\n", LocTypeStr[LocType]);
        fprintf(out, "@    box %g, %g, %g, %g\n", xmin, ymin, xmax, ymax);
        fprintf(out, "@    box linestyle %d\n", LineStyle);
        fprintf(out, "@    box linewidth %d\n", LineWidth);
        fprintf(out, "@    box color %d\n", LineColor);
        fprintf(out, "@    box fill %s\n", BoxFillStr[BoxFill]);
        fprintf(out, "@    box fill color %d\n", BoxColor);
        fprintf(out, "@    box fill pattern %d\n", BoxPattern);
        fprintf(out, "@box def\n");
    }
}

/* reads a line into ptr, adjusting len and renewing ptr if neccesary */
static char *fgets3(FILE *fp, char **ptr, int *len, int maxlen)
{
    int   len_remaining = *len; /* remaining amount of allocated bytes in buf */
    int   curp          = 0;    /* current position in buf to read into */

    do
    {
        if (len_remaining < 2)
        {
            if (*len + STRLEN < maxlen)
            {
                /* This line is longer than len characters, let's increase len! */
                *len          += STRLEN;
                len_remaining += STRLEN;
                srenew(*ptr, *len);
            }
            else
            {
                /*something is wrong, we'll just keep reading and return NULL*/
                len_remaining = STRLEN;
                curp          = 0;
            }
        }
        if (fgets(*ptr + curp, len_remaining, fp) == nullptr)
        {
            /* if last line, skip */
            return nullptr;
        }
        curp         += len_remaining-1; /* overwrite the nul char in next iteration */
        len_remaining = 1;
    }
    while ((std::strchr(*ptr, '\n') == nullptr) && (!feof(fp)));

    if (*len + STRLEN >= maxlen)
    {
        return nullptr; /* this line was too long */
    }

    if (feof(fp))
    {
        /* We reached EOF before '\n', skip this last line. */
        return nullptr;
    }
    {
        /* now remove newline */
        int slen = std::strlen(*ptr);
        if ((*ptr)[slen-1] == '\n')
        {
            (*ptr)[slen-1] = '\0';
        }
    }

    return *ptr;
}

static int wordcount(char *ptr)
{
    int i, n = 0, is[2];
    int cur = 0;
#define prev (1-cur)

    if (nullptr != ptr)
    {
        for (i = 0; (ptr[i] != '\0'); i++)
        {
            is[cur] = std::isspace(ptr[i]);
            if ((0 == i) && !is[cur])
            {
                n++;
            }
            else if ((i > 0)  && (!is[cur] && is[prev]))
            {
                n++;
            }
            cur = prev;
        }
    }
    return n;
}

static char *read_xvgr_string(const char *line)
{
    const char *ptr0, *ptr1;
    char       *str;

    ptr0 = std::strchr(line, '"');
    if (ptr0 != nullptr)
    {
        ptr0++;
        ptr1 = std::strchr(ptr0, '"');
        if (ptr1 != nullptr)
        {
            str            = gmx_strdup(ptr0);
            str[ptr1-ptr0] = '\0';
        }
        else
        {
            str = gmx_strdup("");
        }
    }
    else
    {
        str = gmx_strdup("");
    }

    return str;
}

int read_xvg_legend(const char *fn, double ***y, int *ny,
                    char **subtitle, char ***legend)
{
    FILE    *fp;
    char    *ptr;
    char    *base = nullptr;
    char    *fmt  = nullptr;
    int      k, line = 0, nny, nx, maxx, rval, legend_nalloc, set, nchar;
    double   lf;
    double **yy = nullptr;
    char    *tmpbuf;
    int      len = STRLEN;
    *ny  = 0;
    nny  = 0;
    nx   = 0;
    maxx = 0;
    fp   = gmx_fio_fopen(fn, "r");

    snew(tmpbuf, len);
    if (subtitle != nullptr)
    {
        *subtitle = nullptr;
    }
    legend_nalloc = 0;
    if (legend != nullptr)
    {
        *legend = nullptr;
    }

    while ((ptr = fgets3(fp, &tmpbuf, &len, 10*STRLEN)) != nullptr && ptr[0] != '&')
    {
        line++;
        trim(ptr);
        if (ptr[0] == '@')
        {
            if (legend != nullptr)
            {
                ptr++;
                trim(ptr);
                set = -1;
                if (std::strncmp(ptr, "subtitle", 8) == 0)
                {
                    ptr += 8;
                    if (subtitle != nullptr)
                    {
                        *subtitle = read_xvgr_string(ptr);
                    }
                }
                else if (std::strncmp(ptr, "legend string", 13) == 0)
                {
                    ptr += 13;
                    sscanf(ptr, "%d%n", &set, &nchar);
                    ptr += nchar;
                }
                else if (ptr[0] == 's')
                {
                    ptr++;
                    sscanf(ptr, "%d%n", &set, &nchar);
                    ptr += nchar;
                    trim(ptr);
                    if (std::strncmp(ptr, "legend", 6) == 0)
                    {
                        ptr += 6;
                    }
                    else
                    {
                        set = -1;
                    }
                }
                if (set >= 0)
                {
                    if (set >= legend_nalloc)
                    {
                        legend_nalloc = set + 1;
                        srenew(*legend, legend_nalloc);
                        (*legend)[set] = read_xvgr_string(ptr);
                    }
                }
            }
        }
        else if (ptr[0] != '#')
        {
            if (nny == 0)
            {
                (*ny) = nny = wordcount(ptr);
                /* fprintf(stderr,"There are %d columns in your file\n",nny);*/
                if (nny == 0)
                {
                    return 0;
                }
                snew(yy, nny);
                snew(fmt, 3*nny+1);
                snew(base, 3*nny+1);
            }
            /* Allocate column space */
            if (nx >= maxx)
            {
                maxx += 1024;
                for (k = 0; (k < nny); k++)
                {
                    srenew(yy[k], maxx);
                }
            }
            /* Initiate format string */
            fmt[0]  = '\0';
            base[0] = '\0';

            /* fprintf(stderr,"ptr='%s'\n",ptr);*/
            for (k = 0; (k < nny); k++)
            {
                std::strcpy(fmt, base);
                std::strcat(fmt, "%lf");
                rval = sscanf(ptr, fmt, &lf);
                /* fprintf(stderr,"rval = %d\n",rval);*/
                if ((rval == EOF) || (rval == 0))
                {
                    break;
                }
                yy[k][nx] = lf;
                srenew(fmt, 3*(nny+1)+1);
                srenew(base, 3*nny+1);
                std::strcat(base, "%*s");
            }
            if (k != nny)
            {
                fprintf(stderr, "Only %d columns on line %d in file %s\n",
                        k, line, fn);
                for (; (k < nny); k++)
                {
                    yy[k][nx] = 0.0;
                }
            }
            nx++;
        }
    }
    gmx_fio_fclose(fp);

    *y = yy;
    sfree(tmpbuf);
    sfree(base);
    sfree(fmt);

    if (legend_nalloc > 0)
    {
        if (*ny - 1 > legend_nalloc)
        {
            assert(legend);
            srenew(*legend, *ny-1);
            for (set = legend_nalloc; set < *ny-1; set++)
            {
                (*legend)[set] = nullptr;
            }
        }
    }

    return nx;
}

int read_xvg(const char *fn, double ***y, int *ny)
{
    return read_xvg_legend(fn, y, ny, nullptr, nullptr);
}

void write_xvg(const char *fn, const char *title, int nx, int ny, real **y,
               const char **leg, const gmx_output_env_t *oenv)
{
    FILE *fp;
    int   i, j;

    fp = xvgropen(fn, title, "X", "Y", oenv);
    if (leg)
    {
        xvgr_legend(fp, ny-1, leg, oenv);
    }
    for (i = 0; (i < nx); i++)
    {
        for (j = 0; (j < ny); j++)
        {
            fprintf(fp, "  %12.5e", y[j][i]);
        }
        fprintf(fp, "\n");
    }
    xvgrclose(fp);
}

real **read_xvg_time(const char *fn,
                     gmx_bool bHaveT, gmx_bool bTB, real tb, gmx_bool bTE, real te,
                     int nsets_in, int *nset, int *nval, real *dt, real **t)
{
    FILE      *fp;
#define MAXLINELEN 16384
    char       line0[MAXLINELEN];
    char      *line;
    int        t_nalloc, *val_nalloc, a, narg, n, sin, set, nchar;
    double     dbl;
    gmx_bool   bEndOfSet, bTimeInRange, bFirstLine = TRUE;
    real     **val;

    t_nalloc   = 0;
    *t         = nullptr;
    val        = nullptr;
    val_nalloc = nullptr;
    *nset      = 0;
    *dt        = 0;
    fp         = gmx_fio_fopen(fn, "r");
    for (sin = 0; sin < nsets_in; sin++)
    {
        if (nsets_in == 1)
        {
            narg = 0;
        }
        else
        {
            narg = bHaveT ? 2 : 1;
        }
        n         = 0;
        bEndOfSet = FALSE;
        while (!bEndOfSet && fgets(line0, MAXLINELEN, fp))
        {
            line = line0;
            /* Remove whitespace */
            while (line[0] == ' ' || line[0] == '\t')
            {
                line++;
            }
            bEndOfSet = (line[0] == '&');
            if (line[0] != '#' && line[0] != '@' && line[0] != '\n' && !bEndOfSet)
            {
                if (bFirstLine && bHaveT)
                {
                    /* Check the first line that should contain data */
                    a = sscanf(line, "%lf%lf", &dbl, &dbl);
                    if (a == 0)
                    {
                        gmx_fatal(FARGS, "Expected a number in %s on line:\n%s", fn, line0);
                    }
                    else if (a == 1)
                    {
                        fprintf(stderr, "Found only 1 number on line, "
                                "assuming no time is present.\n");
                        bHaveT = FALSE;
                        if (nsets_in > 1)
                        {
                            narg = 1;
                        }
                    }
                }

                a            = 0;
                bTimeInRange = TRUE;
                while ((a < narg || (nsets_in == 1 && n == 0)) &&
                       sscanf(line, "%lf%n", &dbl, &nchar) == 1 && bTimeInRange)
                {
                    /* Use set=-1 as the time "set" */
                    if (sin)
                    {
                        if (!bHaveT || (a > 0))
                        {
                            set = sin;
                        }
                        else
                        {
                            set = -1;
                        }
                    }
                    else
                    {
                        if (!bHaveT)
                        {
                            set = a;
                        }
                        else
                        {
                            set = a-1;
                        }
                    }
                    if (set == -1 && ((bTB && dbl < tb) || (bTE && dbl > te)))
                    {
                        bTimeInRange = FALSE;
                    }

                    if (bTimeInRange)
                    {
                        if (n == 0)
                        {
                            if (nsets_in == 1)
                            {
                                narg++;
                            }
                            if (set >= 0)
                            {
                                *nset = set+1;
                                srenew(val, *nset);
                                srenew(val_nalloc, *nset);
                                val_nalloc[set] = 0;
                                val[set]        = nullptr;
                            }
                        }
                        if (set == -1)
                        {
                            if (sin == 0)
                            {
                                if (n >= t_nalloc)
                                {
                                    t_nalloc = over_alloc_small(n);
                                    srenew(*t, t_nalloc);
                                }
                                (*t)[n] = dbl;
                            }
                            /* else we should check the time of the next sets with set 0 */
                        }
                        else
                        {
                            if (n >= val_nalloc[set])
                            {
                                val_nalloc[set] = over_alloc_small(n);
                                srenew(val[set], val_nalloc[set]);
                            }
                            val[set][n] = (real)dbl;
                        }
                    }
                    a++;
                    line += nchar;
                }
                if (line0[strlen(line0)-1] != '\n')
                {
                    fprintf(stderr, "File %s does not end with a newline, ignoring the last line\n", fn);
                }
                else if (bTimeInRange)
                {
                    if (a == 0)
                    {
                        fprintf(stderr, "Ignoring invalid line in %s:\n%s", fn, line0);
                    }
                    else
                    {
                        if (a != narg)
                        {
                            fprintf(stderr, "Invalid line in %s:\n%s"
                                    "Using zeros for the last %d sets\n",
                                    fn, line0, narg-a);
                        }
                        n++;
                    }
                }
                if (a > 0)
                {
                    bFirstLine = FALSE;
                }
            }
        }
        if (sin == 0)
        {
            *nval = n;
            if (!bHaveT)
            {
                snew(*t, n);
                for (a = 0; a < n; a++)
                {
                    (*t)[a] = a;
                }
            }
            if (n > 1)
            {
                *dt = (real)((*t)[n-1]-(*t)[0])/(n-1.0);
            }
            else
            {
                *dt = 1;
            }
        }
        else
        {
            if (n < *nval)
            {
                fprintf(stderr, "Set %d is shorter (%d) than the previous set (%d)\n",
                        sin+1, n, *nval);
                *nval = n;
                fprintf(stderr, "Will use only the first %d points of every set\n",
                        *nval);
            }
        }
    }
    gmx_fio_fclose(fp);

    sfree(val_nalloc);

    return val;
}
