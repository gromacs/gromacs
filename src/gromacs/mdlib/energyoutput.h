/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_ENERGYOUTPUT_H
#define GMX_MDLIB_ENERGYOUTPUT_H

#include <cstdio>

#include "gromacs/mdtypes/enerdata.h"

class energyhistory_t;
struct ener_file;
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct SimulationGroups;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct t_ebin;
struct t_expanded;
struct t_fcdata;
struct t_grpopts;
struct t_inputrec;
struct t_lambda;
class t_state;

namespace gmx
{
class Awh;
class Constraints;
}

extern const char *egrp_nm[egNR+1];

/* delta_h block type enum: the kinds of energies written out. */
enum
{
    dhbtDH   = 0, /* delta H BAR energy difference*/
    dhbtDHDL = 1, /* dH/dlambda derivative */
    dhbtEN,       /* System energy */
    dhbtPV,       /* pV term */
    dhbtEXPANDED, /* expanded ensemble statistics */
    dhbtNR
};

namespace gmx
{

// TODO remove use of detail namespace when removing t_mdebin in
// favour of an Impl class.
namespace detail
{
struct t_mdebin;
}

/* The functions & data structures here determine the content for outputting
   the .edr file; the file format and actual writing is done with functions
   defined in enxio.h */

class EnergyOutput
{
    public:
        EnergyOutput();
        /*! \brief Initiate MD energy bin
         *
         * This second phase of construction is needed until we have
         * modules that understand how to request output from
         * EnergyOutput.
         *
         * \todo Refactor to separate a function to write the energy
         * file header. Perhaps transform the remainder into a factory
         * function.
         */
        void prepare(ener_file        *fp_ene,
                     const gmx_mtop_t *mtop,
                     const t_inputrec *ir,
                     FILE             *fp_dhdl,
                     bool              isRerun = false);
        ~EnergyOutput();
        /*! \brief Update the averaging structures.
         *
         * Called every step on which the energies are evaluated. */
        void addDataAtEnergyStep(bool                    bDoDHDL,
                                 bool                    bSum,
                                 double                  time,
                                 real                    tmass,
                                 gmx_enerdata_t         *enerd,
                                 t_state                *state,
                                 t_lambda               *fep,
                                 t_expanded             *expand,
                                 matrix                  lastbox,
                                 tensor                  svir,
                                 tensor                  fvir,
                                 tensor                  vir,
                                 tensor                  pres,
                                 gmx_ekindata_t         *ekind,
                                 rvec                    mu_tot,
                                 const gmx::Constraints *constr);
        /*! \brief Updated the averaging structures
         *
         * Called every step on which the energies are not evaluated.
         *
         * \todo This schedule is known in advance, and should be made
         * an intrinsic behaviour of EnergyOutput, rather than being
         * wastefully called every step. */
        void recordNonEnergyStep();

        /*! \brief Help write quantites to the energy file
         *
         * \todo Perhaps this responsibility should involve some other
         * object visiting all the contributing objects. */
        void printStepToEnergyFile(ener_file *fp_ene, bool bEne, bool bDR, bool bOR,
                                   FILE *log,
                                   int64_t step, double time,
                                   int mode,
                                   t_fcdata *fcd,
                                   SimulationGroups *groups, t_grpopts *opts,
                                   gmx::Awh *awh);
        /*! \brief Get the number of energy terms recorded.
         *
         * \todo Refactor this to return the expected output size,
         * rather than exposing the implementation details about
         * energy terms. */
        int numEnergyTerms() const;
        /*! \brief Getter used for testing t_ebin
         *
         * \todo Find a better approach for this. */
        t_ebin *getEbin();

        /* Between .edr writes, the averages are history dependent,
           and that history needs to be retained in checkpoints.
           These functions set/read the energyhistory_t class
           that is written to checkpoints in checkpoint.c */

        //! Fill the energyhistory_t data.
        void fillEnergyHistory(energyhistory_t * enerhist) const;
        //! Restore from energyhistory_t data.
        void restoreFromEnergyHistory(const energyhistory_t &enerhist);

    private:
        // TODO transform this into an impl class.
        detail::t_mdebin *mdebin = nullptr;
};

} // namespace gmx

//! Open the dhdl file for output
FILE *open_dhdl(const char *filename, const t_inputrec *ir,
                const gmx_output_env_t *oenv);

namespace gmx
{

//! Print an energy-output header to the log file
void print_ebin_header(FILE *log, int64_t steps, double time);

} // namespace gmx

#endif
