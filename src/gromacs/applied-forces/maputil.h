/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMX_APPLIED_FORCES_MAPUTIL_H
#define GMX_APPLIED_FORCES_MAPUTIL_H

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_commrec;
struct t_filenm;
struct t_inputrec;

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/math/griddata/griddata.h"

namespace gmx
{
class Densfit;
struct t_mapdata;

t_mapdata gridDataToMapData(const GridDataReal3D &inputGrid);
GridDataReal3D mapDataToGridData(const t_mapdata &inputMap);


//! \brief Allocate memory for the density data
extern void allocate_density_grid(
        const std::array<int, DIM> &nValues,
        std::vector<float> *vox);

//! \brief If the box turns out to be all 0's, then assign a useful box from the positions
extern void translate_pdb(matrix box, rvec x[], int natoms, gmx_bool bTranslate);

//! \brief Allocate and initialize the atomic weights array TODO
extern void init_weights(gmx_bool bVerbose, Densfit *densfit);


//! \brief Lower-level force calculation routine
extern void do_densfit_forces(Densfit *densfit, matrix box);

//! \brief From the map entries determine and return the grid spacing
extern real get_map_spacing(const t_mapdata &map, FILE *log);

//! \brief Append a number to the output file name
extern void make_filename(const char *outf_name, int ndigit, int file_nr, char newname[STRLEN]);

//! \brief TODO
extern void make_positive(t_mapdata *map_ref);

//! \brief TODO
extern t_mapdata rescale_map(const t_mapdata &map_ref, real scale);

} // namespace gmx

#endif
