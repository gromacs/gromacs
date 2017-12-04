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

#include "densfitmapdata.h"

#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/real.h"

struct t_blocka;
struct t_commrec;
struct t_inputrec;
struct t_filenm;
struct gmx_mtop_t;
struct gmx_domdec_t;
struct gmx_output_env_t;
struct gmx_wallcycle;

namespace gmx
{
class DensfitData;

class LocalDensfitData;

class Densfit
{
    public:
        Densfit();
        Densfit(real sigma, real sigma_dist, real k, /* Spring constant */
                int nat,
                int *ind,                            /* Which of the atoms should be used for spreading */
                real grid_spacing, bool bVerbose);
        Densfit(const DensfitData &parameters);
        ~Densfit();

        /*! \brief Return a copy of the reference map.
         */
        t_mapdata referenceMapCopy() const;

        void makeGroups(char *groupname, t_blocka *grps, char **gnames);

        void setSpacingFromMap(const t_mapdata &map);
        const t_mapdata &simulatedMap();

        //! \brief Compute the forces so that fitting to a cryo-EM density map can be
        //! done
        void do_densfit(real t, gmx_int64_t step, t_inputrec *ir, t_commrec *cr,
                        PaddedArrayRef<RVec> x, matrix box, gmx_wallcycle *wcycle);

        //! \brief Add the density fitting forces to the MD forces and output, return
        //! the correlation coefficient
        real add_forces(rvec *f, gmx_unused t_commrec *cr, gmx_int64_t step,
                        real time);
        //! \brief Calculate the correlation coefficient of the two maps
        real calc_correlation_coeff(FILE *log);

        void setParameters(const DensfitData &parameters);

        void setReferenceMap(const t_mapdata &referenceMap);

        //! \brief Transform the positions x into a density map by replacing each
        //! position by a Gaussian function of width sigma
        void spread_atoms(matrix box);
        //! \brief Make a selection of the home atoms for the density fitting group.
        //! Should be called at every domain decomposition.
        void dd_make_local_df_indices(gmx_domdec_t *dd);
        //! \brief Make a new x array that only contains atoms to be spread
        void assemble_atoms_for_spread(rvec x[]);
        //! \brief Initialize density fitting
        void init(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
                  gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr,
                  const gmx_output_env_t *oenv, gmx_bool bAppend, gmx_bool bVerbose);
        //! \brief Clean up density fitting, close output file
        void finish();
        void broadcast(const t_commrec *cr);
        const DensfitData &parameters() const;
        //! \brief Allocate memory for the simulated density grid
        void new_map_sim_from_ref(const t_mapdata &map_ref);

    private:
        std::unique_ptr<DensfitData> parameters_;
        std::unique_ptr<LocalDensfitData>
        df_; /* Stores non-inputrec density fitting data */

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
        void spread_atoms_low(rvec x[], int natoms, matrix box);
        void do_forces(matrix box);
        FILE *open_out(const char *fn, const gmx_output_env_t *oenv);
        //! \brief Check whether densfit->nstfit is overwritten by environment
        //! variable.
        void get_nstfit_from_env(t_commrec *cr, int nstlist);

        void setPrivateData(rvec *x_densfit_whole, real grid_spacing, bool bVerbose,
                            bool bAppend, bool bParallel);
        void updateParametersThatAreTimeDependent(real time);

        /*! \brief Dumps positions to file.
         * forces are dumped to the default density fitting output file
         */
        void dump_x(int nodeid, gmx_int64_t step);
        /*! \brief Dumps forces to file 'filename'.
         *
         * if no filename given, forces are dumped to the default density fitting
         * output file
         */
        void dump_f(const char *fn, t_commrec *cr);
};

//! \brief Append a number to the output file name
void make_filename(const char *outf_name, int ndigit, int file_nr,
                   char newname[STRLEN]);
}

#endif /* end of include guard: GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFIT_H */
