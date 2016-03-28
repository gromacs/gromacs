/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2017, by the GROMACS development team, led by
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
#ifndef GMX_FILEIO_XVGR_H
#define GMX_FILEIO_XVGR_H

#include <stdio.h>

#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;

/***************************************************
 *            XVGR   DEFINITIONS
 ***************************************************/
enum {
    elNone, elSolid, elDotted, elDashed,
    elLongDashed, elDotDashed, elNR
};
/* xvgr line-styles */

enum {
    ecWhite, ecFrank, ecBlack = ecFrank,
    ecRed, ecGreen, ecBlue, ecYellow, ecBrown, ecGray,
    ecPurple, ecLightBlue, ecViolet, ecHolland, ecLila, ecDarkGray,
    ecAquamarine, ecOlive, ecNR
};
/* xvgr line-colors */

enum {
    eppNone, eppColor, eppPattern, eppNR
};
/* xvgr pattern type */

enum {
    evView, evWorld, evNR
};
/* view type */

/***************************************************
 *            XVGR   ROUTINES
 ***************************************************/

/* Strings such as titles, lables and legends can contain escape sequences
 * for formatting. Currently supported are:
 * \s : start subscript
 * \S : start superscript
 * \N : end sub/superscript
 * \symbol : where symbol is the full name of a greek letter
 *           (see the xvgrstr function in xvgr.c for the full list)
 *           when starting with a capital, a capital symbol will be printed,
 *           note that symbol does not need to be followed by a space
 * \8 : (deprecated) start symbol font
 * \4 : (deprecated) end symbol font
 */

gmx_bool output_env_get_print_xvgr_codes(const struct gmx_output_env_t *oenv);
/* Returns if we should print xmgrace or xmgr codes */

enum {
    exvggtNONE, exvggtXNY, exvggtXYDY, exvggtXYDYDY, exvggtNR
};

void xvgr_header(FILE *fp, const char *title, const std::string &xaxis,
                 const std::string &yaxis, int exvg_graph_type,
                 const struct gmx_output_env_t *oenv);
/* In most cases you want to use xvgropen_type, which does the same thing
 * but takes a filename and opens it.
 */

FILE *xvgropen_type(const char *fn, const char *title, const std::string &xaxis,
                    const std::string &yaxis, int exvg_graph_type,
                    const struct gmx_output_env_t *oenv);
/* Open a file, and write a title, and axis-labels in Xvgr format
 * or write nothing when oenv specifies so.
 * The xvgr graph type enum is defined above.
 */

FILE *xvgropen(const char *fn, const char *title, const std::string &xaxis,
               const std::string &yaxis, const struct gmx_output_env_t *oenv);
/* Calls xvgropen_type with graph type xvggtXNY. */

/* Close xvgr file, and clean up internal file buffers correctly */
void xvgrclose(FILE *fp);

void xvgr_subtitle(FILE *out, const char *subtitle,
                   const struct gmx_output_env_t *oenv);
/* Set the subtitle in xvgr */

void xvgr_view(FILE *out, real xmin, real ymin, real xmax, real ymax,
               const struct gmx_output_env_t *oenv);
/* Set the view in xvgr */

void xvgr_world(FILE *out, real xmin, real ymin, real xmax, real ymax,
                const struct gmx_output_env_t *oenv);
/* Set the world in xvgr */

void xvgrLegend(FILE                           *out,
                const std::vector<std::string> &setNames,
                const struct gmx_output_env_t  *oenv);
/* Make a legend box, and also modifies the view to make room for the legend */

void xvgr_legend(FILE *out, int nsets, const char** setnames,
                 const struct gmx_output_env_t *oenv);
/* Make a legend box, and also modifies the view to make room for the legend */

void xvgr_new_dataset(FILE *out,
                      int nr_first, int nsets, const char **setnames,
                      const struct gmx_output_env_t *oenv);
/* End the previous data set(s) and start new one(s).
    nr_first = the global set number of the first new set (or 0 if no legend)
    nsets = the number of sets (or 0 if no legends)
    setnames = the set names (or NULL if no legends)
 */

void xvgr_line_props(FILE *out, int NrSet, int LineStyle, int LineColor,
                     const struct gmx_output_env_t *oenv);
/* Set xvgr line styles and colors */

void xvgr_box(FILE *out,
              int LocType,
              real xmin, real ymin, real xmax, real ymax,
              int LineStyle, int LineWidth, int LineColor,
              int BoxFill, int BoxColor, int BoxPattern,
              const struct gmx_output_env_t *oenv);
/* Make a box */

int read_xvg_legend(const char *fn, double ***y, int *ny,
                    char **subtitle, char ***legend);
/* Read an xvg file for post processing. The number of rows is returned
 * fn is the filename, y is a pointer to a 2D array (to be allocated by
 * the routine) ny is the number of columns (including X if appropriate).
 * If subtitle!=NULL, read the subtitle (when present),
 * the subtitle string will be NULL when not present.
 * If legend!=NULL, read the legends for the sets (when present),
 * 0 is the first y legend, the legend string will be NULL when not present.
 */

int read_xvg(const char *fn, double ***y, int *ny);
/* As read_xvg_legend, but does not read legends. */

void write_xvg(const char *fn, const char *title, int nx, int ny, real **y,
               const char** leg, const struct gmx_output_env_t *oenv);
/* Write a two D array (y) of dimensions nx rows times
 * ny columns to a file. If leg != NULL it will be written too.
 */


/* This function reads ascii (xvg) files and extracts the data sets to a
 * two dimensional array which is returned.
 */
real **read_xvg_time(const char *fn,
                     gmx_bool bHaveT,
                     gmx_bool bTB, real tb,
                     gmx_bool bTE, real te,
                     int nsets_in, int *nset, int *nval,
                     real *dt, real **t);
#endif
