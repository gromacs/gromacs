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
struct t_densfit;
struct t_filenm;
struct t_inputrec;
struct t_mapdata;

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/cstringutil.h"

namespace gmx
{

//! \brief Make a structure containing density fitting data
extern t_densfit new_t_densfit(
        real      sigma,
        real      sigma_dist,
        real      k,        /* Spring constant */
        int       nat,
        int      *ind,      /* Which of the atoms should be used for spreading */
        real      grid_spacing,
        gmx_bool  bVerbose);

//! \brief Print the header information stored in file on disk
extern void map_print_header(FILE *fpout, const char *fn);

//! \brief Dumps positions to file with file name "fn", if no filename given, forces are dumped to the default density fitting output file
extern void dump_x(t_densfit *densfit, t_commrec *cr, gmx_int64_t step);

//! \brief Dumps forces to file with file name "fn", if no filename given, forces are dumped to the default density fitting output file
extern void dump_f(const char *fn, t_densfit *densfit, t_commrec *cr);

//! \brief Allocate memory for the density data
extern void allocate_density_grid(
        const int           grid[3],
        std::vector<float> *vox);

//! \brief Allocate memory for the simulated density grid
extern void new_map_sim_from_ref(
        t_densfit       *densfit,
        const t_mapdata &map_ref);

//! \brief If the box turns out to be all 0's, then assign a useful box from the positions
extern void translate_pdb(matrix box, rvec x[], int natoms, gmx_bool bTranslate);

//! \brief If bRead==TRUE : Write map data to fn; FALSE: Read  map data from fn and store in map, allocate memory
extern void do_map_ccp4(gmx_bool bRead, t_mapdata *map, const char *fn,
                        gmx_bool bOverwrite, gmx_bool bVerbose, FILE *DumpHeader);

//! \brief Calculate the correlation coefficient of the two maps
extern real calc_correlation_coeff(t_densfit *densfit, FILE *log);

//! \brief Transform the positions x into a density map by replacing each position by a Gaussian function of width sigma
extern void spread_atoms(t_densfit  *densfit, matrix box);

//! \brief Initialize density fitting
extern void init_density_fitting(
        FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
        gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv,
        gmx_bool bAppendFiles, gmx_bool bVerbose);

//! \brief Clean up density fitting, close output file
extern void finish_density_fitting(const t_densfit &densfit);

//! \brief Allocate and initialize the atomic weights array TODO
extern void init_weights(gmx_bool bVerbose, t_densfit *densfit);

//! \brief Make a new x array that only contains atoms to be spread
extern void assemble_atoms_for_spread(
        t_densfit *densfit,
        rvec       x[]);

//! \brief Make a selection of the home atoms for the density fitting group. Should be called at every domain decomposition.
extern void dd_make_local_df_indices(gmx_domdec_t *dd, t_densfit *densfit);

//! \brief Compute the forces so that fitting to a cryo-EM density map can be done
extern void do_densfit(
        real                           t,
        gmx_int64_t                    step,
        gmx_bool                       bOutputMap,
        t_inputrec                    *ir,
        t_commrec                     *cr,
        gmx::PaddedArrayRef<gmx::RVec> x,
        matrix                         box,
        gmx_wallcycle                 *wcycle);

//! \brief Lower-level force calculation routine
extern void do_densfit_forces(t_densfit *densfit, matrix box);

//! \brief Add the density fitting forces to the MD forces and output, return the correlation coefficient
extern real add_densfit_forces(t_inputrec *ir, rvec *f, gmx_unused t_commrec *cr, gmx_int64_t step, real time);

//! \brief From the map entries determine and return the grid spacing
extern real get_map_spacing(const t_mapdata &map, FILE *log);

//! \brief TODO
extern void couple_map_spacing_to_box(matrix box, t_densfit *densfit);

//! \brief Append a number to the output file name
extern void make_filename(const char *outf_name, int ndigit, int file_nr, char newname[STRLEN]);

//! \brief TODO
extern void make_positive(t_mapdata *map_ref);

//! \brief TODO
extern t_mapdata rescale_map(const t_mapdata &map_ref, real scale);

} // namespace gmx

#endif
