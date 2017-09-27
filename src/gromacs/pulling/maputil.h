/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
#ifndef _maputil_h
#define _maputil_h

struct gmx_domdec_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_filenm;

#include "gromacs/utility/cstringutil.h"

extern t_densfit *new_t_densfit(
        real      sigma,
        real      sigma_dist,
        real      k,        /* Spring constant */
        int       nat,
        int      *ind,      /* Which of the atoms should be used for spreading */
        real      grid_spacing,
        gmx_bool  bVerbose);

extern void gmx_map_print_header(FILE *fpout, const char *fn);

extern void dump_x(t_densfit *densfit, t_commrec *cr, gmx_int64_t step);
extern void dump_f(const char *fn, t_densfit *densfit, t_commrec *cr);

extern void allocate_density_grid(
        int grid[3],
        float **vox);

/* Allocate memory for the simulated density grid */
extern void new_map_sim_from_ref(
        t_densfit *densfit,
        t_mapdata *map_ref);

extern void translate_pdb(matrix box, rvec x[], int natoms, gmx_bool bTranslate);

extern void gmx_do_map_ccp4(gmx_bool bRead, t_mapdata **ptr_map, const char *fn,
                            gmx_bool bOverwrite, gmx_bool bVerbose, FILE *DumpHeader);

/* Calculate the correlation coefficient of the two maps */
extern real calc_correlation_coeff(t_densfit *densfit, FILE *log);

/* Transform the positions x into a density map by replacing
 * each position by a Gaussian function of width sigma */
extern void spread_atoms(t_densfit  *densfit, matrix box);

extern void do_densfit_forces(t_densfit *densfit, matrix box);

/* Initialize density fitting */
extern void init_density_fitting(
        FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
        gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv,
        gmx_bool bAppendFiles, gmx_bool bVerbose);

/* Allocate and initialize the atomic weights array */
extern void init_weights(gmx_bool bVerbose, t_densfit *densfit);

/* If Intels vector math library is present, set the precision here */
extern void gmx_set_vml_precision(FILE *fp);

/* Make a new x array that only contains atoms to be spread */
extern void assemble_atoms_for_spread(
        t_densfit *densfit,
        rvec       x[]);

extern void dd_make_local_df_indices(gmx_domdec_t *dd, t_densfit *densfit);
/* Make a selection of the home atoms for the density fitting group.
 * Should be called at every domain decomposition. */

/* Compute the forces so that fitting to a cryo-EM density map can be done */
extern void do_densfit(
        real           t,
        gmx_int64_t    step,
        gmx_bool       bOutputMap,
        t_inputrec    *ir,
        t_commrec     *cr,
        rvec           x[],
        matrix         box,
        gmx_wallcycle *wcycle);

extern void add_densfit_forces(t_inputrec *ir, rvec *f, t_commrec *cr, gmx_int64_t step, real time);

/* From the map entries determine and return the grid spacing */
extern real get_map_spacing(t_mapdata *map, FILE *log);

extern void couple_map_spacing_to_box(matrix box, t_densfit *densfit);

/* Append a number to the output file name */
extern void make_filename(const char *outf_name, int ndigit, int file_nr, char newname[STRLEN]);

extern void make_positive(t_mapdata *map_ref);

extern t_mapdata *rescale_map(t_mapdata *map_ref, real scale);

#endif
