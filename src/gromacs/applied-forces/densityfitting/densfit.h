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

#ifndef GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFIT_H
#define GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFIT_H

#include <array>
#include <vector>
#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/paddedvector.h"

struct t_blocka;
struct t_inpfile;
struct t_commrec;
struct t_inputrec;
struct t_filenm;
struct gmx_mtop_t;
struct gmx_domdec_t;
struct gmx_output_env_t;
struct gmx_wallcycle;

namespace gmx
{

struct gmxDensfit;

/* This data structure holds all the data extracted from a ccp4 cmap datafile */
struct t_mapdata
{
    char                    *title;
    int                      datamode;
    std::array<float, 2*DIM> cell;
    std::array<int, DIM>     grid;
    std::array<int, DIM>     origin;
    std::array<int, DIM>     axes_order;
    std::array<int, DIM>     map_dim;
    int                      spacegroup = 1;
    float                    min        = 0.;
    float                    max        = 0.;
    double                   mean       = 0.;
    double                   rms        = 0.;
    std::array<float, 9>     skew_mat   = {};
    std::array<float, DIM>   skew_trans = {};
    std::vector<float>       vox; /* 1d continuous array for 3d voxel density        */
};

/* Abstract type for density fitting mdrun-data only defined in maputil.cpp    */

class Densfit
{
    public:
        Densfit();
        Densfit(real      sigma,
                real      sigma_dist,
                real      k,    /* Spring constant */
                int       nat,
                int      *ind,  /* Which of the atoms should be used for spreading */
                real      grid_spacing,
                bool      bVerbose);
        ~Densfit();
        void setPrivateData(rvec      *x_densfit_whole, /* Can be NULL!                       */
                            real       grid_spacing,
                            bool       bVerbose,
                            bool       bAppend,
                            bool       bParallel);
        void makeGroups(char *groupname, t_blocka *grps, char **gnames);
        void updateParametersThatAreTimeDependent(real time);
        //! \brief Dumps positions to file with file name "fn", if no filename given, forces are dumped to the default density fitting output file
        void dump_x(int nodeid, gmx_int64_t step);
        //! \brief Dumps forces to file with file name "fn", if no filename given, forces are dumped to the default density fitting output file
        void dump_f(const char *fn, t_commrec * cr);
        void couple_map_spacing_to_box(matrix box);
        int   npoints;                                  /* Number of time points for the refinement with
                                                           different values for sigma, k and T             */
        real *time_values;                              /* Array (size npoints) of time values (ps)        */
        real *temp_values;                              /* Same for T (K)                                  */
        real *k_values;                                 /* Same for k (kJ/mol)                             */
        real *sigma_values;                             /* Same for sigma (nm)                             */
        real  dist;                                     /* Cutoff distance (in sigmas) for density spreading
                                                           spread each atom only in the range
                                                           -dist ... +dist, ignore contributions further
                                                           away                                            */
        int                         nstfit;             /* Only recalculate V_fit every nstfit time steps  */
        int                         nstout;             /* Write diagnostic output every nstout time steps */
        int                         nstmapout;          /* Keep simulated maps in these intervals          */
        int                         nat;                /* Number of atoms to be spread on the map         */
        int                        *ind;                /* The global atoms numbers                        */
        int                         nat_loc;            /* Number of local spreading atoms                 */
        int                         nalloc_loc;         /* Allocation size for ind_loc and weight_loc      */
        int                        *ind_loc;            /* Local spreading indices                         */
        int                         nweight;            /* The number of weights (0 or nat)                */
        real                       *weight;             /* Weights (use all 1 when weight==NULL)           */
        real                       *weight_loc;         /* Weights for the local indices                   */

        t_mapdata                   map_ref;            /* The reference=experimental map to fit to        */
        t_mapdata                   map_sim;            /* The map simulated from atomic position data     */
        gmx_bool                    bKeepAndNumberMaps; /* ... or just keep the last one ...               */
        std::unique_ptr<gmxDensfit> df;                 /* Stores non-inputrec density fitting data        */

        void do_densfit(real t, gmx_int64_t step, gmx_bool bOutputMap, t_inputrec *ir, t_commrec *cr, PaddedArrayRef< RVec> x, matrix box, gmx_wallcycle *wcycle);

        real add_forces( rvec *f, gmx_unused t_commrec *cr, gmx_int64_t step, real time);
        //! \brief Calculate the correlation coefficient of the two maps
        real calc_correlation_coeff(FILE *log);
        //! \brief Transform the positions x into a density map by replacing each position by a Gaussian function of width sigma
        void spread_atoms(matrix box);
        //! \brief Make a selection of the home atoms for the density fitting group. Should be called at every domain decomposition.
        void dd_make_local_df_indices(gmx_domdec_t *dd);
        //! \brief Make a new x array that only contains atoms to be spread
        void assemble_atoms_for_spread(rvec       x[]);
        //! \brief Initialize density fitting
        void init(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
                  gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv,
                  gmx_bool bAppend, gmx_bool bVerbose);
        //! \brief Clean up density fitting, close output file
        void finish();

    private:
        real interpolateLinearly(real time, int npoints, real* time_values, real* y_values);
        /*! \brief Calculate a simulated density
         *
         * Make the density map, assuming that the protein is whole, and in the first
         * quadrant. Transform the positions x into a density map by replacing
         * each position by a Gaussian function of width sigma
         *
         * \param[in] x  Atom positions
         * \param[in] natoms number of Atom positions
         * \param[in] box the simulation box
         */
        void spread_atoms_low( rvec x[], int natoms, matrix box);
        void do_forces(matrix box);
        FILE *open_out(const char *fn, const gmx_output_env_t *oenv);
        //! \brief Check whether densfit->nstfit is overwritten by environment variable.
        void get_nstfit_from_env(t_commrec *cr, int nstlist);



};

}

#endif /* end of include guard: GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFIT_H */
