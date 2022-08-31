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
#ifndef GMX_FILEIO_XVGR_H
#define GMX_FILEIO_XVGR_H

#include <stdio.h>

#include <optional>
#include <string>
#include <vector>

#include "gromacs/math/multidimarray.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;
namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx
/***************************************************
 *            XVGR   DEFINITIONS
 ***************************************************/
enum
{
    elNone,
    elSolid,
    elDotted,
    elDashed,
    elLongDashed,
    elDotDashed,
    elNR
};
/* xvgr line-styles */

enum
{
    ecWhite,
    ecFrank,
    ecBlack = ecFrank,
    ecRed,
    ecGreen,
    ecBlue,
    ecYellow,
    ecBrown,
    ecGray,
    ecPurple,
    ecLightBlue,
    ecViolet,
    ecHolland,
    ecLila,
    ecDarkGray,
    ecAquamarine,
    ecOlive,
    ecNR
};
/* xvgr line-colors */

enum
{
    eppNone,
    eppColor,
    eppPattern,
    eppNR
};
/* xvgr pattern type */

enum
{
    evView,
    evWorld,
    evNR
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

gmx_bool output_env_get_print_xvgr_codes(const struct gmx_output_env_t* oenv);
/* Returns if we should print xmgrace or xmgr codes */

enum
{
    exvggtNONE,
    exvggtXNY,
    exvggtXYDY,
    exvggtXYDYDY,
    exvggtNR
};

void xvgr_header(FILE*                          fp,
                 const char*                    title,
                 const std::string&             xaxis,
                 const std::string&             yaxis,
                 int                            exvg_graph_type,
                 const struct gmx_output_env_t* oenv);
/* In most cases you want to use xvgropen_type, which does the same thing
 * but takes a filename and opens it.
 */

FILE* xvgropen_type(const char*                    fn,
                    const char*                    title,
                    const std::string&             xaxis,
                    const std::string&             yaxis,
                    int                            exvg_graph_type,
                    const struct gmx_output_env_t* oenv);
/* Open a file, and write a title, and axis-labels in Xvgr format
 * or write nothing when oenv specifies so.
 * The xvgr graph type enum is defined above.
 */

FILE* xvgropen(const char*                    fn,
               const char*                    title,
               const std::string&             xaxis,
               const std::string&             yaxis,
               const struct gmx_output_env_t* oenv);
/* Calls xvgropen_type with graph type xvggtXNY. */

/* Close xvgr file, and clean up internal file buffers correctly */
void xvgrclose(FILE* fp);

void xvgr_subtitle(FILE* out, const char* subtitle, const struct gmx_output_env_t* oenv);
/* Set the subtitle in xvgr */

void xvgr_view(FILE* out, real xmin, real ymin, real xmax, real ymax, const struct gmx_output_env_t* oenv);
/* Set the view in xvgr */

void xvgr_world(FILE* out, real xmin, real ymin, real xmax, real ymax, const struct gmx_output_env_t* oenv);
/* Set the world in xvgr */

//! Prepare a legend box, also modifies the view to make room for the legend
void xvgrLegend(FILE* out, gmx::ArrayRef<const std::string> setNames, const struct gmx_output_env_t* oenv);

/*! \brief
 * End the previous data set(s) and start new one(s).
 *
 * \param[in] out File to write to.
 * \param[in] firstSetNumber Global number of the first data set, or 0 if no legend.
 * \param[in] setNames View on collection of strings for legend in the data set.
 * \param[in] oenv Global output enivornment handling.
 */
void xvgrNewDataset(FILE*                            out,
                    int                              firstSetNumber,
                    gmx::ArrayRef<const std::string> setNames,
                    const struct gmx_output_env_t*   oenv);

void xvgr_line_props(FILE* out, int NrSet, int LineStyle, int LineColor, const struct gmx_output_env_t* oenv);
/* Set xvgr line styles and colors */

void xvgr_box(FILE*                          out,
              int                            LocType,
              real                           xmin,
              real                           ymin,
              real                           xmax,
              real                           ymax,
              int                            LineStyle,
              int                            LineWidth,
              int                            LineColor,
              int                            BoxFill,
              int                            BoxColor,
              int                            BoxPattern,
              const struct gmx_output_env_t* oenv);
/* Make a box */

int read_xvg_legend(const char* fn, double*** y, int* ny, char** subtitle, char*** legend);
/* Read an xvg file for post processing. The number of rows is returned
 * fn is the filename, y is a pointer to a 2D array (to be allocated by
 * the routine) ny is the number of columns (including X if appropriate).
 * If subtitle!=NULL, read the subtitle (when present),
 * the subtitle string will be NULL when not present.
 * If legend!=NULL, read the legends for the sets (when present),
 * 0 is the first y legend, the legend string will be NULL when not present.
 */

/* \brief Read only the data from an xvg file for post processing.
 *
 * Note: this function is deprecated in favor of readXvg, which is
 *       used under the hood in this function.
 *
 * \param[out]    nx Number of rows.
 * \param[in]     fn Xvg file to read.
 * \param[in/out] y  Pointer to 2D array (allocated by the routine).
 * \param[in/out] ny Number of columns.
 *
 * Todo: Port all read_xvg calls to use readXvgData
 */
int read_xvg(const char* fn, double*** y, int* ny);

/* \brief Read only the data from an xvg file for post processing.
 *
 * \param[out] XvgData Data in row major.
 * \param[in]  fn      Xvg file to read.
 */
gmx::MultiDimArray<std::vector<double>, gmx::dynamicExtents2D> readXvgData(const std::string& fn);


void write_xvg(const char*                      fn,
               const char*                      title,
               int                              nx,
               int                              ny,
               real**                           y,
               gmx::ArrayRef<const std::string> leg,
               const struct gmx_output_env_t*   oenv);
/* Write a two D array (y) of dimensions nx rows times
 * ny columns to a file. If leg != NULL it will be written too.
 */

/*! \brief
 * Read xvg data as a time series.
 *
 * Allows truncation of data series to exclude time points.
 * Expects first row to be
 * time series. Only one set can be read in at the same time.
 *
 * \returns Data series in row major, first row being the time series.
 * \param[in] fn Xvg file to read
 * \param[in] startTime Optional first time to read.
 * \param[in] endTime Optional last time to read.
 */
gmx::MultiDimArray<std::vector<double>, gmx::dynamicExtents2D>
readXvgTimeSeries(const std::string& fn, std::optional<real> startTime, std::optional<real> endTime);


/*!\brief
 *  This function reads ascii (xvg) files and extracts the data sets to a
 * two dimensional array which is returned.
 *
 * NOTE: This function is deprecated and shouldn't be used for new code.
 */
real** read_xvg_time(const char* fn,
                     gmx_bool    bHaveT,
                     gmx_bool    bTB,
                     real        tb,
                     gmx_bool    bTE,
                     real        te,
                     int         nsets_in,
                     int*        nset,
                     int*        nval,
                     real*       dt,
                     real**      t);
#endif
