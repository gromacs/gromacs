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

#include "matio.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <optional>
#include <regex>
#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx;

static const char mapper[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/"
        "?";
#define NMAP static_cast<long int>(sizeof(mapper) / sizeof(mapper[0]))

real** mk_matrix(int nx, int ny, gmx_bool b1D)
{
    int    i;
    real** m;

    snew(m, nx);
    if (b1D)
    {
        snew(m[0], nx * ny);
    }

    for (i = 0; (i < nx); i++)
    {
        if (b1D)
        {
            m[i] = &(m[0][i * ny]);
        }
        else
        {
            snew(m[i], ny);
        }
    }

    return m;
}

void done_matrix(int nx, real*** m)
{
    int i;

    for (i = 0; (i < nx); i++)
    {
        sfree((*m)[i]);
    }
    sfree(*m);
    *m = nullptr;
}

t_matelmt searchcmap(ArrayRef<const t_mapping> map, t_xpmelmt c)
{
    auto findIt = std::find_if(map.begin(), map.end(), [&c](const auto& m) { return (m.code == c); });
    return (findIt == map.end()) ? -1 : (findIt - map.begin());
}

//! Read the mapping table from in, return number of entries
static std::vector<t_mapping> getcmap(FILE* in, const std::filesystem::path& fn)
{
    int                    i, n;
    char                   line[STRLEN];
    char                   code[STRLEN], desc[STRLEN];
    double                 r, g, b;
    std::vector<t_mapping> m;

    if (fgets2(line, STRLEN - 1, in) == nullptr)
    {
        gmx_fatal(FARGS,
                  "Not enough lines in colormap file %s"
                  "(just wanted to read number of entries)",
                  fn.string().c_str());
    }
    sscanf(line, "%d", &n);
    m.resize(n);
    for (i = 0; (i < n); i++)
    {
        if (fgets2(line, STRLEN - 1, in) == nullptr)
        {
            gmx_fatal(FARGS,
                      "Not enough lines in colormap file %s"
                      "(should be %d, found only %d)",
                      fn.string().c_str(),
                      n + 1,
                      i);
        }
        sscanf(line, "%s%s%lf%lf%lf", code, desc, &r, &g, &b);
        m[i].code.c1 = code[0];
        m[i].code.c2 = 0;
        m[i].desc    = desc;
        m[i].rgb.r   = r;
        m[i].rgb.g   = g;
        m[i].rgb.b   = b;
    }

    return m;
}

std::vector<t_mapping> readcmap(const std::filesystem::path& fn)
{
    FilePtr in = openLibraryFile(fn);
    return getcmap(in.get(), fn);
}

void printcmap(FILE* out, int n, const t_mapping map[])
{
    int i;

    fprintf(out, "%d\n", n);
    for (i = 0; (i < n); i++)
    {
        fprintf(out,
                "%c%c  %20s  %10g  %10g  %10g\n",
                map[i].code.c1 ? map[i].code.c1 : ' ',
                map[i].code.c2 ? map[i].code.c2 : ' ',
                map[i].desc.c_str(),
                map[i].rgb.r,
                map[i].rgb.g,
                map[i].rgb.b);
    }
}

void writecmap(const std::filesystem::path& fn, int n, const t_mapping map[])
{
    FILE* out;

    out = gmx_fio_fopen(fn, "w");
    printcmap(out, n, map);
    gmx_fio_fclose(out);
}

static char* fgetline(char** line, int llmax, int* llalloc, FILE* in)
{
    char* fg;

    if (llmax > *llalloc)
    {
        srenew(*line, llmax + 1);
        *llalloc = llmax;
    }
    fg = fgets(*line, llmax, in);
    trim(*line);

    return fg;
}

static void skipstr(char* line)
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
        line[c - i] = line[c];
        c++;
    }
    line[c - i] = '\0';
}

static char* line2string(char** line)
{
    int i;

    if (*line != nullptr)
    {
        while (((*line)[0] != '\"') && ((*line)[0] != '\0'))
        {
            (*line)++;
        }

        if ((*line)[0] != '\"')
        {
            return nullptr;
        }
        (*line)++;

        i = 0;
        while (((*line)[i] != '\"') && ((*line)[i] != '\0'))
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

//! If a label named \c label is found in \c line, return it. Otherwise return empty string.
static std::optional<std::string> findLabelInLine(const std::string& line, const std::string& label)
{
    std::regex  re(".*\\s" + label + ":[\\s]*\"(.*)\"");
    std::smatch match;
    if (std::regex_search(line, match, re) && match.size() > 1)
    {
        return std::make_optional<std::string>(match.str(1));
    }
    return std::nullopt;
}

//! Read and return a matrix from \c in
static t_matrix read_xpm_entry(FILE* in)
{
    char *                     line_buf = nullptr, *line = nullptr, *str;
    std::optional<std::string> title, legend, xLabel, yLabel, matrixType;
    int                        i, m, col_len, nch = 0, llmax;
    int                        llalloc = 0;
    unsigned int               r, g, b;
    double                     u;
    gmx_bool                   bGetOnWithIt, bSetLine;
    t_xpmelmt                  c;

    llmax = STRLEN;

    while ((nullptr != fgetline(&line_buf, llmax, &llalloc, in))
           && (std::strncmp(line_buf, "static", 6) != 0))
    {
        std::string lineString = line_buf;
        if (!title.has_value())
        {
            title = findLabelInLine(lineString, "title");
        }
        if (!legend.has_value())
        {
            legend = findLabelInLine(lineString, "legend");
        }
        if (!xLabel.has_value())
        {
            xLabel = findLabelInLine(lineString, "x-label");
        }
        if (!yLabel.has_value())
        {
            yLabel = findLabelInLine(lineString, "y-label");
        }
        if (!matrixType.has_value())
        {
            matrixType = findLabelInLine(lineString, "type");
        }
    }

    if (!line_buf || strncmp(line_buf, "static", 6) != 0)
    {
        gmx_input("Invalid XPixMap");
    }

    t_matrix mm;

    mm.title   = title.value_or("");
    mm.legend  = legend.value_or("");
    mm.label_x = xLabel.value_or("");
    mm.label_y = yLabel.value_or("");

    if (matrixType.has_value() && (gmx_strcasecmp(matrixType->c_str(), "Discrete") == 0))
    {
        mm.bDiscrete = TRUE;
    }

    if (debug)
    {
        fprintf(debug,
                "%s %s %s %s\n",
                mm.title.c_str(),
                mm.legend.c_str(),
                mm.label_x.c_str(),
                mm.label_y.c_str());
    }

    /* Read sizes */
    int nmap     = 0;
    bGetOnWithIt = FALSE;
    while (!bGetOnWithIt && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)))
    {
        line = line_buf;
        while ((line[0] != '\"') && (line[0] != '\0'))
        {
            line++;
        }

        if (line[0] == '\"')
        {
            line2string(&line);
            sscanf(line, "%d %d %d %d", &(mm.nx), &(mm.ny), &nmap, &nch);
            if (nch > 2)
            {
                gmx_fatal(FARGS, "Sorry can only read xpm's with at most 2 characters per pixel\n");
            }
            if (mm.nx <= 0 || mm.ny <= 0)
            {
                gmx_fatal(FARGS, "Dimensions of xpm-file have to be larger than 0\n");
            }
            llmax        = std::max(STRLEN, mm.nx + 10);
            bGetOnWithIt = TRUE;
        }
    }
    if (debug)
    {
        fprintf(debug, "mm.nx %d mm.ny %d nmap %d nch %d\n", mm.nx, mm.ny, nmap, nch);
    }
    if (nch == 0)
    {
        gmx_fatal(FARGS, "Number of characters per pixel not found in xpm\n");
    }

    /* Read color map */
    mm.map.resize(nmap);
    m = 0;
    while ((m < nmap) && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)))
    {
        line = std::strchr(line_buf, '\"');
        if (line)
        {
            line++;
            /* Read xpm color map entry */
            mm.map[m].code.c1 = line[0];
            if (nch == 1)
            {
                mm.map[m].code.c2 = 0;
            }
            else
            {
                mm.map[m].code.c2 = line[1];
            }
            line += nch;
            str = std::strchr(line, '#');
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
                    mm.map[m].rgb.r = r / 255.0;
                    mm.map[m].rgb.g = g / 255.0;
                    mm.map[m].rgb.b = b / 255.0;
                }
                else if (col_len == 12)
                {
                    sscanf(line, "%*s #%4x%4x%4x", &r, &g, &b);
                    mm.map[m].rgb.r = r / 65535.0;
                    mm.map[m].rgb.g = g / 65535.0;
                    mm.map[m].rgb.b = b / 65535.0;
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
                mm.map[m].rgb.r = 1;
                mm.map[m].rgb.g = 1;
                mm.map[m].rgb.b = 1;
            }
            line = std::strchr(line, '\"');
            line++;
            line2string(&line);
            mm.map[m].desc = line;
            m++;
        }
    }
    if (m != nmap)
    {
        gmx_fatal(FARGS,
                  "Number of read colors map entries (%d) does not match the number in the header "
                  "(%d)",
                  m,
                  nmap);
    }

    /* Read axes, if there are any */
    bSetLine = FALSE;
    do
    {
        if (bSetLine)
        {
            line = line_buf;
        }
        bSetLine = TRUE;
        GMX_RELEASE_ASSERT(line, "Need to have valid line to parse");
        if (strstr(line, "x-axis"))
        {
            line = std::strstr(line, "x-axis");
            skipstr(line);
            mm.axis_x.reserve(mm.nx + 1);
            while (sscanf(line, "%lf", &u) == 1)
            {
                if (gmx::ssize(mm.axis_x) > mm.nx)
                {
                    gmx_fatal(FARGS, "Too many x-axis labels in xpm (max %d)", mm.nx);
                }
                else if (gmx::ssize(mm.axis_x) == mm.nx)
                {
                    mm.flags |= MAT_SPATIAL_X;
                }
                mm.axis_x.push_back(u);
                skipstr(line);
            }
        }
        else if (std::strstr(line, "y-axis"))
        {
            line = std::strstr(line, "y-axis");
            skipstr(line);
            mm.axis_y.reserve(mm.ny + 1);
            while (sscanf(line, "%lf", &u) == 1)
            {
                if (gmx::ssize(mm.axis_y) > mm.ny)
                {
                    gmx_fatal(FARGS, "Too many y-axis labels in xpm (max %d)", mm.ny);
                }
                else if (gmx::ssize(mm.axis_y) == mm.ny)
                {
                    mm.flags |= MAT_SPATIAL_Y;
                }
                mm.axis_y.push_back(u);
                skipstr(line);
            }
        }
    } while ((line[0] != '\"') && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)));

    /* Read matrix */
    mm.matrix.resize(mm.nx, mm.ny);
    int rowIndex = mm.ny - 1;
    bSetLine     = FALSE;
    do
    {
        if (bSetLine)
        {
            line = line_buf;
        }
        bSetLine = TRUE;
        if (rowIndex % (1 + mm.ny / 100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100 * (mm.ny - rowIndex)) / mm.ny);
        }
        while ((line[0] != '\"') && (line[0] != '\0'))
        {
            line++;
        }
        if (line[0] != '\"')
        {
            gmx_fatal(FARGS, "Not enough characters in row %d of the matrix\n", rowIndex + 1);
        }
        else
        {
            line++;
            for (i = 0; i < mm.nx; i++)
            {
                c.c1 = line[nch * i];
                if (nch == 1)
                {
                    c.c2 = 0;
                }
                else
                {
                    c.c2 = line[nch * i + 1];
                }
                mm.matrix(i, rowIndex) = searchcmap(mm.map, c);
            }
            rowIndex--;
        }
    } while ((rowIndex >= 0) && (nullptr != fgetline(&line_buf, llmax, &llalloc, in)));
    if (rowIndex >= 0)
    {
        gmx_incons("Not enough rows in the matrix");
    }

    sfree(line_buf);
    return mm;
}

std::vector<t_matrix> read_xpm_matrix(const std::filesystem::path& fnm)
{
    FILE* in;
    char* line    = nullptr;
    int   llalloc = 0;

    in = gmx_fio_fopen(fnm, "r");

    std::vector<t_matrix> mat;
    while (nullptr != fgetline(&line, STRLEN, &llalloc, in))
    {
        if (std::strstr(line, "/* XPM */"))
        {
            mat.emplace_back(read_xpm_entry(in));
        }
    }
    gmx_fio_fclose(in);

    if (mat.empty())
    {
        gmx_file("Invalid XPixMap");
    }

    sfree(line);

    return mat;
}

real** matrix2real(const t_matrix* in, real** out)
{
    double tmp;

    std::vector<real> rmap(in->map.size());

    for (gmx::Index i = 0; i != gmx::ssize(in->map); ++i)
    {
        if ((in->map[i].desc.empty()) || (sscanf(in->map[i].desc.c_str(), "%lf", &tmp) != 1))
        {
            fprintf(stderr,
                    "Could not convert matrix to reals,\n"
                    "color map entry %zd has a non-real description: \"%s\"\n",
                    i,
                    in->map[i].desc.c_str());
            return nullptr;
        }
        rmap[i] = tmp;
    }

    if (out == nullptr)
    {
        snew(out, in->nx);
        for (int i = 0; i < in->nx; i++)
        {
            snew(out[i], in->ny);
        }
    }
    for (int i = 0; i < in->nx; i++)
    {
        for (int j = 0; j < in->ny; j++)
        {
            out[i][j] = rmap[in->matrix(i, j)];
        }
    }

    fprintf(stderr, "Converted a %dx%d matrix with %zu levels to reals\n", in->nx, in->ny, in->map.size());

    return out;
}

static void write_xpm_header(FILE*              out,
                             const std::string& title,
                             const std::string& legend,
                             const std::string& label_x,
                             const std::string& label_y,
                             gmx_bool           bDiscrete)
{
    fprintf(out, "/* XPM */\n");
    fprintf(out, "/* This file can be converted to EPS by the GROMACS program xpm2ps */\n");
    fprintf(out, "/* title:   \"%s\" */\n", title.c_str());
    fprintf(out, "/* legend:  \"%s\" */\n", legend.c_str());
    fprintf(out, "/* x-label: \"%s\" */\n", label_x.c_str());
    fprintf(out, "/* y-label: \"%s\" */\n", label_y.c_str());
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
    return std::min(std::max(0, static_cast<int>(((mid - lo) / (hi - lo)) * (nlevels - 1))), nlevels - 1);
}

static void
write_xpm_map3(FILE* out, int n_x, int n_y, int* nlevels, real lo, real mid, real hi, t_rgb rlo, t_rgb rmid, t_rgb rhi)
{
    int    i, nmid;
    double r, g, b, clev_lo, clev_hi;

    if (*nlevels > NMAP * NMAP)
    {
        fprintf(stderr,
                "Warning, too many levels (%d) in matrix, using %d only\n",
                *nlevels,
                static_cast<int>(NMAP * NMAP));
        *nlevels = NMAP * NMAP;
    }
    else if (*nlevels < 2)
    {
        fprintf(stderr, "Warning, too few levels (%d) in matrix, using 2 instead\n", *nlevels);
        *nlevels = 2;
    }
    if (!((mid >= lo) && (mid < hi)))
    {
        gmx_fatal(FARGS, "Lo: %f, Mid: %f, Hi: %f\n", lo, mid, hi);
    }

    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %d %d\",\n", n_x, n_y, *nlevels, (*nlevels <= NMAP) ? 1 : 2);

    nmid    = calc_nmid(*nlevels, lo, mid, hi);
    clev_lo = nmid;
    clev_hi = (*nlevels - 1 - nmid);
    for (i = 0; (i < nmid); i++)
    {
        r = rlo.r + (i * (rmid.r - rlo.r) / clev_lo);
        g = rlo.g + (i * (rmid.g - rlo.g) / clev_lo);
        b = rlo.b + (i * (rmid.b - rlo.b) / clev_lo);
        fprintf(out,
                "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[i % NMAP],
                (*nlevels <= NMAP) ? ' ' : mapper[i / NMAP],
                static_cast<unsigned int>(std::round(255 * r)),
                static_cast<unsigned int>(std::round(255 * g)),
                static_cast<unsigned int>(std::round(255 * b)),
                ((nmid - i) * lo + i * mid) / clev_lo);
    }
    for (i = 0; (i < (*nlevels - nmid)); i++)
    {
        r = rmid.r + (i * (rhi.r - rmid.r) / clev_hi);
        g = rmid.g + (i * (rhi.g - rmid.g) / clev_hi);
        b = rmid.b + (i * (rhi.b - rmid.b) / clev_hi);
        fprintf(out,
                "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[(i + nmid) % NMAP],
                (*nlevels <= NMAP) ? ' ' : mapper[(i + nmid) / NMAP],
                static_cast<unsigned int>(std::round(255 * r)),
                static_cast<unsigned int>(std::round(255 * g)),
                static_cast<unsigned int>(std::round(255 * b)),
                ((*nlevels - 1 - nmid - i) * mid + i * hi) / clev_hi);
    }
}

static void pr_simple_cmap(FILE* out, real lo, real hi, int nlevel, t_rgb rlo, t_rgb rhi, int i0)
{
    int  i;
    real r, g, b, fac;

    for (i = 0; (i < nlevel); i++)
    {
        fac = (i + 1.0) / (nlevel);
        r   = rlo.r + fac * (rhi.r - rlo.r);
        g   = rlo.g + fac * (rhi.g - rlo.g);
        b   = rlo.b + fac * (rhi.b - rlo.b);
        fprintf(out,
                "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[(i + i0) % NMAP],
                (nlevel <= NMAP) ? ' ' : mapper[(i + i0) / NMAP],
                static_cast<unsigned int>(std::round(255 * r)),
                static_cast<unsigned int>(std::round(255 * g)),
                static_cast<unsigned int>(std::round(255 * b)),
                lo + fac * (hi - lo));
    }
}

static void pr_discrete_cmap(FILE* out, int* nlevel, int i0)
{
    t_rgb rgbd[16] = {
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

    int i, n;

    *nlevel = std::min(16, *nlevel);
    n       = *nlevel;
    for (i = 0; (i < n); i++)
    {
        fprintf(out,
                "\"%c%c c #%02X%02X%02X \" /* \"%3d\" */,\n",
                mapper[(i + i0) % NMAP],
                (n <= NMAP) ? ' ' : mapper[(i + i0) / NMAP],
                static_cast<unsigned int>(std::round(255 * rgbd[i].r)),
                static_cast<unsigned int>(std::round(255 * rgbd[i].g)),
                static_cast<unsigned int>(std::round(255 * rgbd[i].b)),
                i);
    }
}


static void write_xpm_map_split(FILE*      out,
                                int        n_x,
                                int        n_y,
                                const int* nlevel_top,
                                real       lo_top,
                                real       hi_top,
                                t_rgb      rlo_top,
                                t_rgb      rhi_top,
                                gmx_bool   bDiscreteColor,
                                int*       nlevel_bot,
                                real       lo_bot,
                                real       hi_bot,
                                t_rgb      rlo_bot,
                                t_rgb      rhi_bot)
{
    int ntot;

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


static void write_xpm_map(FILE* out, int n_x, int n_y, int* nlevels, real lo, real hi, t_rgb rlo, t_rgb rhi)
{
    int  i, nlo;
    real invlevel, r, g, b;

    if (*nlevels > NMAP * NMAP)
    {
        fprintf(stderr,
                "Warning, too many levels (%d) in matrix, using %d only\n",
                *nlevels,
                static_cast<int>(NMAP * NMAP));
        *nlevels = NMAP * NMAP;
    }
    else if (*nlevels < 2)
    {
        fprintf(stderr, "Warning, too few levels (%d) in matrix, using 2 instead\n", *nlevels);
        *nlevels = 2;
    }

    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %d %d\",\n", n_x, n_y, *nlevels, (*nlevels <= NMAP) ? 1 : 2);

    invlevel = 1.0 / (*nlevels - 1);
    for (i = 0; (i < *nlevels); i++)
    {
        nlo = *nlevels - 1 - i;
        r   = (nlo * rlo.r + i * rhi.r) * invlevel;
        g   = (nlo * rlo.g + i * rhi.g) * invlevel;
        b   = (nlo * rlo.b + i * rhi.b) * invlevel;
        fprintf(out,
                "\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                mapper[i % NMAP],
                (*nlevels <= NMAP) ? ' ' : mapper[i / NMAP],
                static_cast<unsigned int>(std::round(255 * r)),
                static_cast<unsigned int>(std::round(255 * g)),
                static_cast<unsigned int>(std::round(255 * b)),
                (nlo * lo + i * hi) * invlevel);
    }
}

static void writeXpmAxis(FILE* out, const char* axis, ArrayRef<const real> label)
{
    if (label.empty())
    {
        return;
    }
    for (gmx::Index i = 0; i != ssize(label); ++i)
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

static void write_xpm_data(FILE* out, int n_x, int n_y, const real* const* mat, real lo, real hi, int nlevels)
{
    int  i, j, c;
    real invlevel;

    invlevel = (nlevels - 1) / (hi - lo);
    for (j = n_y - 1; (j >= 0); j--)
    {
        if (j % (1 + n_y / 100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100 * (n_y - j)) / n_y);
        }
        fprintf(out, "\"");
        for (i = 0; (i < n_x); i++)
        {
            c = roundToInt((mat[i][j] - lo) * invlevel);
            if (c < 0)
            {
                c = 0;
            }
            if (c >= nlevels)
            {
                c = nlevels - 1;
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

static void write_xpm_data3(FILE* out, int n_x, int n_y, real** mat, real lo, real mid, real hi, int nlevels)
{
    int  i, j, c = 0, nmid;
    real invlev_lo, invlev_hi;

    nmid      = calc_nmid(nlevels, lo, mid, hi);
    invlev_hi = (nlevels - 1 - nmid) / (hi - mid);
    invlev_lo = (nmid) / (mid - lo);

    for (j = n_y - 1; (j >= 0); j--)
    {
        if (j % (1 + n_y / 100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100 * (n_y - j)) / n_y);
        }
        fprintf(out, "\"");
        for (i = 0; (i < n_x); i++)
        {
            if (mat[i][j] >= mid)
            {
                c = nmid + roundToInt((mat[i][j] - mid) * invlev_hi);
            }
            else if (mat[i][j] >= lo)
            {
                c = roundToInt((mat[i][j] - lo) * invlev_lo);
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
                c = nlevels - 1;
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

static void write_xpm_data_split(FILE*  out,
                                 int    n_x,
                                 int    n_y,
                                 real** mat,
                                 real   lo_top,
                                 real   hi_top,
                                 int    nlevel_top,
                                 real   lo_bot,
                                 real   hi_bot,
                                 int    nlevel_bot)
{
    int  i, j, c;
    real invlev_top, invlev_bot;

    invlev_top = (nlevel_top - 1) / (hi_top - lo_top);
    invlev_bot = (nlevel_bot - 1) / (hi_bot - lo_bot);

    for (j = n_y - 1; (j >= 0); j--)
    {
        if (j % (1 + n_y / 100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100 * (n_y - j)) / n_y);
        }
        fprintf(out, "\"");
        for (i = 0; (i < n_x); i++)
        {
            if (i < j)
            {
                c = nlevel_bot + roundToInt((mat[i][j] - lo_top) * invlev_top);
                if ((c < nlevel_bot) || (c >= nlevel_bot + nlevel_top))
                {
                    gmx_fatal(FARGS,
                              "Range checking i = %d, j = %d, c = %d, bot = %d, top = %d "
                              "matrix[i,j] = %f",
                              i,
                              j,
                              c,
                              nlevel_bot,
                              nlevel_top,
                              mat[i][j]);
                }
            }
            else if (i > j)
            {
                c = roundToInt((mat[i][j] - lo_bot) * invlev_bot);
                if ((c < 0) || (c >= nlevel_bot + nlevel_bot))
                {
                    gmx_fatal(FARGS,
                              "Range checking i = %d, j = %d, c = %d, bot = %d, top = %d "
                              "matrix[i,j] = %f",
                              i,
                              j,
                              c,
                              nlevel_bot,
                              nlevel_top,
                              mat[i][j]);
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

void write_xpm_m(FILE* out, t_matrix m)
{
    gmx_bool  bOneChar;
    t_xpmelmt c;

    bOneChar = (m.map[0].code.c2 == 0);
    write_xpm_header(out, m.title, m.legend, m.label_x, m.label_y, m.bDiscrete);
    fprintf(out, "static char *gromacs_xpm[] = {\n");
    fprintf(out, "\"%d %d   %zu %d\",\n", m.nx, m.ny, m.map.size(), bOneChar ? 1 : 2);
    for (const auto& map : m.map)
    {
        fprintf(out,
                "\"%c%c c #%02X%02X%02X \" /* \"%s\" */,\n",
                map.code.c1,
                bOneChar ? ' ' : map.code.c2,
                static_cast<unsigned int>(std::round(map.rgb.r * 255)),
                static_cast<unsigned int>(std::round(map.rgb.g * 255)),
                static_cast<unsigned int>(std::round(map.rgb.b * 255)),
                map.desc.c_str());
    }
    writeXpmAxis(out, "x", m.axis_x);
    writeXpmAxis(out, "y", m.axis_y);
    for (int j = m.ny - 1; (j >= 0); j--)
    {
        if (j % (1 + m.ny / 100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (100 * (m.ny - j)) / m.ny);
        }
        fprintf(out, "\"");
        if (bOneChar)
        {
            for (int i = 0; (i < m.nx); i++)
            {
                fprintf(out, "%c", m.map[m.matrix(i, j)].code.c1);
            }
        }
        else
        {
            for (int i = 0; (i < m.nx); i++)
            {
                c = m.map[m.matrix(i, j)].code;
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

void write_xpm3(FILE*              out,
                unsigned int       flags,
                const std::string& title,
                const std::string& legend,
                const std::string& label_x,
                const std::string& label_y,
                int                n_x,
                int                n_y,
                real               axis_x[],
                real               axis_y[],
                real*              mat[],
                real               lo,
                real               mid,
                real               hi,
                t_rgb              rlo,
                t_rgb              rmid,
                t_rgb              rhi,
                int*               nlevels)
{
    /* See write_xpm.
     * Writes a colormap varying as rlo -> rmid -> rhi.
     */

    if (hi <= lo)
    {
        gmx_fatal(FARGS, "hi (%g) <= lo (%g)", hi, lo);
    }

    write_xpm_header(out, title, legend, label_x, label_y, FALSE);
    write_xpm_map3(out, n_x, n_y, nlevels, lo, mid, hi, rlo, rmid, rhi);
    writeXpmAxis(out, "x", ArrayRef<real>(axis_x, axis_x + n_x + ((flags & MAT_SPATIAL_X) != 0U ? 1 : 0)));
    writeXpmAxis(out, "y", ArrayRef<real>(axis_y, axis_y + n_y + ((flags & MAT_SPATIAL_Y) != 0U ? 1 : 0)));
    write_xpm_data3(out, n_x, n_y, mat, lo, mid, hi, *nlevels);
}

void write_xpm3(FILE*                                          out,
                unsigned int                                   flags,
                const std::string&                             title,
                const std::string&                             legend,
                const std::string&                             label_x,
                const std::string&                             label_y,
                gmx::ArrayRef<real>                            axis_x,
                gmx::ArrayRef<real>                            axis_y,
                gmx::basic_mdspan<real, gmx::dynamicExtents2D> mat,
                real                                           lo,
                real                                           mid,
                real                                           hi,
                t_rgb                                          rlo,
                t_rgb                                          rmid,
                t_rgb                                          rhi,
                int*                                           nlevels)
{
    real** tempMatrix;
    snew(tempMatrix, mat.extent(0));
    for (int i = 0; i < mat.extent(0); ++i)
    {
        snew(tempMatrix[i], mat.extent(1));
        for (int j = 0; j < mat.extent(1); ++j)
        {
            tempMatrix[i][j] = mat(i, j);
        }
    }
    write_xpm3(out,
               flags,
               title,
               legend,
               label_x,
               label_y,
               axis_x.size(),
               axis_y.size(),
               axis_x.data(),
               axis_y.data(),
               tempMatrix,
               lo,
               mid,
               hi,
               rlo,
               rmid,
               rhi,
               nlevels);
    for (int i = 0; i < mat.extent(0); ++i)
    {
        sfree(tempMatrix[i]);
    }
    sfree(tempMatrix);
}


void write_xpm_split(FILE*              out,
                     unsigned int       flags,
                     const std::string& title,
                     const std::string& legend,
                     const std::string& label_x,
                     const std::string& label_y,
                     int                n_x,
                     int                n_y,
                     real               axis_x[],
                     real               axis_y[],
                     real*              mat[],
                     real               lo_top,
                     real               hi_top,
                     int*               nlevel_top,
                     t_rgb              rlo_top,
                     t_rgb              rhi_top,
                     real               lo_bot,
                     real               hi_bot,
                     int*               nlevel_bot,
                     gmx_bool           bDiscreteColor,
                     t_rgb              rlo_bot,
                     t_rgb              rhi_bot)
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

    write_xpm_header(out, title, legend, label_x, label_y, FALSE);
    write_xpm_map_split(
            out, n_x, n_y, nlevel_top, lo_top, hi_top, rlo_top, rhi_top, bDiscreteColor, nlevel_bot, lo_bot, hi_bot, rlo_bot, rhi_bot);
    writeXpmAxis(out, "x", ArrayRef<real>(axis_x, axis_x + n_x + ((flags & MAT_SPATIAL_X) != 0U ? 1 : 0)));
    writeXpmAxis(out, "y", ArrayRef<real>(axis_y, axis_y + n_y + ((flags & MAT_SPATIAL_Y) != 0U ? 1 : 0)));
    write_xpm_data_split(out, n_x, n_y, mat, lo_top, hi_top, *nlevel_top, lo_bot, hi_bot, *nlevel_bot);
}

void write_xpm_split(FILE*                                          out,
                     unsigned int                                   flags,
                     const std::string&                             title,
                     const std::string&                             legend,
                     const std::string&                             label_x,
                     const std::string&                             label_y,
                     gmx::ArrayRef<real>                            axis_x,
                     gmx::ArrayRef<real>                            axis_y,
                     gmx::basic_mdspan<real, gmx::dynamicExtents2D> mat,
                     real                                           lo_top,
                     real                                           hi_top,
                     int*                                           nlevel_top,
                     t_rgb                                          rlo_top,
                     t_rgb                                          rhi_top,
                     real                                           lo_bot,
                     real                                           hi_bot,
                     int*                                           nlevel_bot,
                     gmx_bool                                       bDiscreteColor,
                     t_rgb                                          rlo_bot,
                     t_rgb                                          rhi_bot)
{
    real** tempMatrix;
    snew(tempMatrix, mat.extent(0));
    for (int i = 0; i < mat.extent(0); ++i)
    {
        snew(tempMatrix[i], mat.extent(1));
        for (int j = 0; j < mat.extent(1); ++j)
        {
            tempMatrix[i][j] = mat(i, j);
        }
    }
    write_xpm_split(out,
                    flags,
                    title,
                    legend,
                    label_x,
                    label_y,
                    axis_x.size(),
                    axis_y.size(),
                    axis_x.data(),
                    axis_y.data(),
                    tempMatrix,
                    lo_top,
                    hi_top,
                    nlevel_top,
                    rlo_top,
                    rhi_top,
                    lo_bot,
                    hi_bot,
                    nlevel_bot,
                    bDiscreteColor,
                    rlo_bot,
                    rhi_bot);
    for (int i = 0; i < mat.extent(0); ++i)
    {
        sfree(tempMatrix[i]);
    }
    sfree(tempMatrix);
}


void write_xpm(FILE*              out,
               unsigned int       flags,
               const std::string& title,
               const std::string& legend,
               const std::string& label_x,
               const std::string& label_y,
               int                n_x,
               int                n_y,
               const real         axis_x[],
               const real         axis_y[],
               const real* const  mat[],
               real               lo,
               real               hi,
               t_rgb              rlo,
               t_rgb              rhi,
               int*               nlevels)
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

    write_xpm_header(out, title, legend, label_x, label_y, FALSE);
    write_xpm_map(out, n_x, n_y, nlevels, lo, hi, rlo, rhi);
    writeXpmAxis(out,
                 "x",
                 ArrayRef<const real>(axis_x, axis_x + n_x + ((flags & MAT_SPATIAL_X) != 0U ? 1 : 0)));
    writeXpmAxis(out,
                 "y",
                 ArrayRef<const real>(axis_y, axis_y + n_y + ((flags & MAT_SPATIAL_Y) != 0U ? 1 : 0)));
    write_xpm_data(out, n_x, n_y, mat, lo, hi, *nlevels);
}

void write_xpm(FILE*                                                out,
               unsigned int                                         flags,
               const std::string&                                   title,
               const std::string&                                   legend,
               const std::string&                                   label_x,
               const std::string&                                   label_y,
               gmx::ArrayRef<const real>                            axis_x,
               gmx::ArrayRef<const real>                            axis_y,
               gmx::basic_mdspan<const real, gmx::dynamicExtents2D> mat,
               real                                                 lo,
               real                                                 hi,
               t_rgb                                                rlo,
               t_rgb                                                rhi,
               int*                                                 nlevels)
{
    real** tempMatrix;
    snew(tempMatrix, mat.extent(0));
    for (int i = 0; i < mat.extent(0); ++i)
    {
        snew(tempMatrix[i], mat.extent(1));
        for (int j = 0; j < mat.extent(1); ++j)
        {
            tempMatrix[i][j] = mat(i, j);
        }
    }
    write_xpm(out,
              flags,
              title,
              legend,
              label_x,
              label_y,
              axis_x.size(),
              axis_y.size(),
              axis_x.data(),
              axis_y.data(),
              tempMatrix,
              lo,
              hi,
              rlo,
              rhi,
              nlevels);
    for (int i = 0; i < mat.extent(0); ++i)
    {
        sfree(tempMatrix[i]);
    }
    sfree(tempMatrix);
}
