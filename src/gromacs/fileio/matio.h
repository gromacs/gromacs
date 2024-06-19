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

#ifndef GMX_FILEIO_MATIO_H
#define GMX_FILEIO_MATIO_H

#include <cstdio>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/fileio/rgb.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx

/*! \brief Models an XPM element
 *
 * \libinternal */
struct t_xpmelmt
{
    //! Should all be non-zero (and printable and not '"')
    char c1 = 0;
    /*! \brief Should all be zero (single char color names: smaller xpm's)
     * or should all be non-zero (double char color names: more colors). */
    char c2 = 0;
};

//! Type of a matrix element
typedef short t_matelmt;

/*! \brief Maps an XPM element to an RGB color and a string description.
 *
 * \libinternal */
struct t_mapping
{
    //! XPM element code
    t_xpmelmt code;
    //! Description
    std::string desc;
    //! RGB color
    t_rgb rgb;
};

//! Comparison for mapping entry.
inline bool operator==(const t_xpmelmt& e1, const t_xpmelmt& e2)
{
    return (e1.c1 == e2.c1) && (e1.c2 == e2.c2);
}

#define MAT_SPATIAL_X (1 << 0)
#define MAT_SPATIAL_Y (1 << 1)
/* Defines if x and y are spatial dimensions,
 * when not, there are n axis ticks at the middle of the elements,
 * when set, there are n+1 axis ticks at the edges of the elements.
 */

//! Convenience typedef
template<typename T>
using DynamicMatrix2D =
        gmx::MultiDimArray<std::vector<T>, gmx::extents<gmx::dynamic_extent, gmx::dynamic_extent>>;

/*! \brief A matrix of integers, plus supporting values, such as used in XPM output
 *
 * \libinternal
 *
 * \todo nx and ny should probably be replaced by operations on the
 * extent of matrix, but currently there is only limited ability to
 * iterate over contents of matrix. */
struct t_matrix
{
    //! Defines if x and y are spatial dimensions. See comments on MAT_*.
    unsigned int flags = 0;
    //! Size of the matrix in x
    int nx = 0;
    //! Size of the matrix in y
    int ny = 0;
    //! Matrix title
    std::string title;
    //! Label for the continuous legend
    std::string legend;
    //! Label for the x-axis
    std::string label_x;
    //! Label for the y-axis
    std::string label_y;
    //! Whether the quantity described is discrete or continuous.
    bool bDiscrete = false;
    //! The x-ticklabels
    std::vector<real> axis_x;
    //! The y-ticklabels
    std::vector<real> axis_y;
    //! Element x,y is matrix(x, y)
    DynamicMatrix2D<t_matelmt> matrix;
    //! Color levels for the output(?)
    std::vector<t_mapping> map;
};

//! Seach in the \c map for code \c c and return its optional entry, or -1 if not found.
t_matelmt searchcmap(gmx::ArrayRef<const t_mapping> map, t_xpmelmt c);

//! Read the mapping table from fn, return number of entries
std::vector<t_mapping> readcmap(const std::filesystem::path& fn);

void printcmap(FILE* out, int n, const t_mapping map[]);
/* print mapping table to out */

void writecmap(const std::filesystem::path& fn, int n, const t_mapping map[]);
/* print mapping table to fn */

//! Reads and returns a number of matrices from .xpm file \c fnm.
std::vector<t_matrix> read_xpm_matrix(const std::filesystem::path& fnm);

real** matrix2real(const t_matrix* in, real** out);
/* Converts an matrix in a t_matrix struct to a matrix of reals
 * When mat==NULL memory will be allocated
 * Returns NULL when something went wrong
 */

void write_xpm_m(FILE* out, t_matrix m);
/* Writes a t_matrix struct to .xpm file */

/*! \brief
 * See write_xpm.
 *
 * Writes a colormap varying as rlo -> rmid -> rhi.
 */
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
                int*               nlevels);

//! \copydoc write_xpm3
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
                int*                                           nlevels);
/*! \brief
 * See write_xpm.
 *
 * Writes a colormap with separate above and below diagonal colormaps.
 * If \c bDiscrete then a colormap with 16 fixed colors is used, first of
 * which is white.
 */
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
                     t_rgb              rhi_bot);

//! \copydoc write_xpm_split
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
                     t_rgb                                          rhi_bot);


/*! \brief Writes matrix data to XPM file.
 * \param[in] out        xpm file
 * \param[in] flags      flags, defined types/matrix.h
 *                       MAT_SPATIAL_X
 *                       MAT_SPATIAL_Y
 *                       Defines if x and y are spatial dimensions,
 *                       when not, there are n axis ticks at the middle of the elements,
 *                       when set, there are n+1 axis ticks at the edges of the elements.
 * \param[in] title      matrix title
 * \param[in] legend     label for the continuous legend
 * \param[in] label_x    label for the x-axis
 * \param[in] label_y    label for the y-axis
 * \param[in] n_x, n_y   size of the matrix
 * \param[in] axis_x[]   the x-ticklabels (n_x or n_x+1)
 * \param[in] axis_y[]   the y-ticklables (n_y or n_y+1)
 * \param[in] *mat[]     element x,y is mat[x][y]
 * \param[in] lo         output lower than lo is set to lo
 * \param[in] hi         output higher than hi is set to hi
 * \param[in] rlo        rgb value for level lo
 * \param[in] rhi        rgb value for level hi
 * \param[in] nlevels    number of color levels for the output
 */
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
               int*               nlevels);

//! \copydoc write_xpm
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
               int*                                                 nlevels);

real** mk_matrix(int nx, int ny, gmx_bool b1D);

void done_matrix(int nx, real*** m);

#endif /* GMX_FILEIO_MATIO_H */
