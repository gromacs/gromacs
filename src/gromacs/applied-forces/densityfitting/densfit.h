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
class LocalAtomSet;
class LocalAtomSetManager;
class LocalDensfitData;

class DensfitOutput
{
    public:
        DensfitOutput(FILE * fplog, const char *fnLogFile, const char * fnDensity, const gmx_output_env_t *oenv, bool bAppend);
        ~DensfitOutput();
        FILE * logFile();
        std::string outputMapFileName() const;
        std::string outputMapFileNameWithStep(gmx_int64_t step) const;
    private:
        bool        bAppend_;
        FILE      * logFile_;
        std::string fnDensity_;        /**< Output filename for simulated density maps     */
};


//! \brief Append a number to the output file name
void make_filename(const char *outf_name, int ndigit, int file_nr,
                   char newname[STRLEN]);
}


struct Densfit
{
    public:
        Densfit();
        Densfit(const gmx::DensfitData &parameters);
        ~Densfit();

        void initOutputFromCommandLineParameters(FILE *fplog, int nfile, const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bAppend, gmx_bool bVerbose);
        /*! \brief Return a copy of the reference map.
         */
        gmx::GridDataReal3D        referenceMapCopy() const;

        const gmx::GridDataReal3D &simulatedMap();

        //! \brief Compute the forces so that fitting to a cryo-EM density map can be
        //! done
        void do_densfit(real t, gmx_int64_t step, const t_commrec *cr,
                        gmx::PaddedArrayRef<gmx::RVec> x, matrix box, gmx_wallcycle *wcycle);

        //! \brief Add the density fitting forces to the MD forces and output.
        void add_forces(rvec *f, gmx_int64_t step,
                        real time);

        /*! \brief Indicate that the group's shift vectors for this structure need to be
         * updated at the next call to communicate_group_positions.
         */
        void                    triggerShiftUpdate();
        //! \brief Initialize density fitting
        void                    init(int ePBC, gmx_mtop_t *mtop, rvec *x, matrix box, gmx::LocalAtomSetManager * atomsets, const t_commrec *cr);
        void                    broadcast(const t_commrec *cr);
        const gmx::DensfitData &parameters() const;

    private:
        std::unique_ptr<gmx::DensfitData>      parameters_; //< constant, global densityfitting data
        std::unique_ptr<gmx::LocalDensfitData> df_;         //< variable, node-local data
        std::unique_ptr<gmx::DensfitOutput>    output_;     //< parameters for outputting data
        static const std::string               name_;

        void                                   get_nstfit_from_env(const t_commrec *cr, int nstlist);

        void                                   setLocalData(rvec *x_densfit_whole, real grid_spacing, bool bVerbose,
                                                            bool bAppend, bool bParallel);

};

#endif /* end of include guard: GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFIT_H */
