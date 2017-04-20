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

#ifndef GMX_FILEIO_MATIO_H
#define GMX_FILEIO_MATIO_H

#include <stdio.h>

#include <string>

#include "gromacs/fileio/rgb.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

typedef struct {
    char c1; /* should all be non-zero (and printable and not '"') */
    char c2; /*
              * should all be zero (single char color names: smaller xpm's)
              * or should all be non-zero (double char color names: more colors)
              */
} t_xpmelmt;

typedef short t_matelmt;

typedef struct {
    t_xpmelmt   code; /* see comment for t_xpmelmt */
    const char *desc;
    t_rgb       rgb;
} t_mapping;

#define MAT_SPATIAL_X (1<<0)
#define MAT_SPATIAL_Y (1<<1)
/* Defines if x and y are spatial dimensions,
 * when not, there are n axis ticks at the middle of the elements,
 * when set, there are n+1 axis ticks at the edges of the elements.
 */

typedef struct {
    unsigned int flags; /* The possible flags are defined above */
    int          nx, ny;
    int          y0;
    char         title[256];
    char         legend[256];
    char         label_x[256];
    char         label_y[256];
    gmx_bool     bDiscrete;
    real        *axis_x;
    real        *axis_y;
    t_matelmt  **matrix;
    int          nmap;
    t_mapping   *map;
} t_matrix;
/* title      matrix title
 * legend     label for the continuous legend
 * label_x    label for the x-axis
 * label_y    label for the y-axis
 * nx, ny     size of the matrix
 * axis_x[]   the x-ticklabels
 * axis_y[]   the y-ticklables
 * *matrix[]  element x,y is matrix[x][y]
 * nmap       number of color levels for the output(?)
 */

gmx_bool matelmt_cmp(t_xpmelmt e1, t_xpmelmt e2);

t_matelmt searchcmap(int n, t_mapping map[], t_xpmelmt c);
/* Seach in the map for code 'c' and return entry number.
 * return -1 if not found
 */

int getcmap(FILE *in, const char *fn, t_mapping **map);
/* Read the mapping table from in, return number of entries */

int readcmap(const char *fn, t_mapping **map);
/* Read the mapping table from fn, return number of entries */

void printcmap(FILE *out, int n, t_mapping map[]);
/* print mapping table to out */

void writecmap(const char *fn, int n, t_mapping map[]);
/* print mapping table to fn */

int read_xpm_matrix(const char *fnm, t_matrix **mat);
/* Reads a number of matrices from .xpm file fnm and returns this number */

real **matrix2real(t_matrix *in, real **out);
/* Converts an matrix in a t_matrix struct to a matrix of reals
 * When mat==NULL memory will be allocated
 * Returns NULL when something went wrong
 */

void write_xpm_m(FILE *out, t_matrix m);
/* Writes a t_matrix struct to .xpm file */

void write_xpm3(FILE *out, unsigned int flags,
                const std::string &title, const std::string &legend,
                const std::string &label_x, const std::string &label_y,
                int n_x, int n_y, real axis_x[], real axis_y[],
                real *mat[], real lo, real mid, real hi,
                t_rgb rlo, t_rgb rmid, t_rgb rhi, int *nlevels);
/* See write_xpm.
 * Writes a colormap varying as rlo -> rmid -> rhi.
 */
void write_xpm_split(FILE *out, unsigned int flags,
                     const std::string &title, const std::string &legend,
                     const std::string &label_x, const std::string &label_y,
                     int n_x, int n_y, real axis_x[], real axis_y[],
                     real *mat[],
                     real lo_top, real hi_top, int *nlevel_top,
                     t_rgb rlo_top, t_rgb rhi_top,
                     real lo_bot, real hi_bot, int *nlevel_bot,
                     gmx_bool bDiscreteColor,
                     t_rgb rlo_bot, t_rgb rhi_bot);
/* See write_xpm.
 * Writes a colormap with separate above and below diagonal colormaps.
 * If bDiscrete then a colormap with 16 fixed colors is used, first  of
 * which is white.
 */

void write_xpm(FILE *out, unsigned int flags,
               const std::string &title, const std::string &legend,
               const std::string &label_x, const std::string &label_y,
               int n_x, int n_y, real t_x[], real t_y[],
               real *mat[], real lo, real hi,
               t_rgb rlo, t_rgb rhi, int *nlevels);
/* out        xpm file
 * flags      flags, defined types/matrix.h
 *            MAT_SPATIAL_X
 *            MAT_SPATIAL_Y
 *            Defines if x and y are spatial dimensions,
 *            when not, there are n axis ticks at the middle of the elements,
 *            when set, there are n+1 axis ticks at the edges of the elements.
 * title      matrix title
 * legend     label for the continuous legend
 * label_x    label for the x-axis
 * label_y    label for the y-axis
 * n_x, n_y   size of the matrix
 * axis_x[]   the x-ticklabels (n_x or n_x+1)
 * axis_y[]   the y-ticklables (n_y or n_y+1)
 * *mat[]     element x,y is mat[x][y]
 * lo         output lower than lo is set to lo
 * hi         output higher than hi is set to hi
 * rlo        rgb value for level lo
 * rhi        rgb value for level hi
 * nlevels    number of color levels for the output
 */

real **mk_matrix(int nx, int ny, gmx_bool b1D);

void done_matrix(int nx, real ***m);

#endif  /* GMX_FILEIO_MATIO_H */
