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
#include "gmxpre.h"

#include "gromacs/pulling/output.h"

#include <cstdio>

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/mdtypes/pullhistory.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

static std::string append_before_extension(const std::string& pathname, const std::string& to_append)
{
    /* Appends to_append before last '.' in pathname */
    size_t extPos = pathname.find_last_of('.');
    if (extPos == std::string::npos)
    {
        return pathname + to_append;
    }
    else
    {
        return pathname.substr(0, extPos) + to_append + pathname.substr(extPos, std::string::npos);
    }
}

static void addToPullxHistory(pull_t* pull)
{
    GMX_RELEASE_ASSERT(pull->coordForceHistory, "Pull history does not exist.");
    pull->coordForceHistory->numValuesInXSum++;
    for (size_t c = 0; c < pull->coord.size(); c++)
    {
        const pull_coord_work_t& pcrd        = pull->coord[c];
        PullCoordinateHistory&   pcrdHistory = pull->coordForceHistory->pullCoordinateSums[c];

        pcrdHistory.value += pcrd.spatialData.value;
        pcrdHistory.valueRef += pcrd.value_ref;

        for (int m = 0; m < DIM; m++)
        {
            pcrdHistory.dr01[m] += pcrd.spatialData.dr01[m];
            pcrdHistory.dr23[m] += pcrd.spatialData.dr23[m];
            pcrdHistory.dr45[m] += pcrd.spatialData.dr45[m];
        }
        if (pcrd.params_.eGeom == PullGroupGeometry::Cylinder)
        {
            for (int m = 0; m < DIM; m++)
            {
                pcrdHistory.dynaX[m] += pcrd.dynamicGroup0->x[m];
            }
        }
    }
    for (size_t g = 0; g < pull->group.size(); g++)
    {
        PullGroupHistory& pgrpHistory = pull->coordForceHistory->pullGroupSums[g];
        for (int m = 0; m < DIM; m++)
        {
            pgrpHistory.x[m] += pull->group[g].x[m];
        }
    }
}

static void addToPullfHistory(pull_t* pull)
{
    GMX_RELEASE_ASSERT(pull->coordForceHistory, "Pull history does not exist.");
    pull->coordForceHistory->numValuesInFSum++;
    for (size_t c = 0; c < pull->coord.size(); c++)
    {
        const pull_coord_work_t& pcrd = pull->coord[c];

        PullCoordinateHistory& pcrdHistory = pull->coordForceHistory->pullCoordinateSums[c];

        pcrdHistory.scalarForce += pcrd.scalarForce;
    }
}

static void pullResetHistory(PullHistory* history, bool resetXHistory, bool resetFHistory)
{
    if (resetXHistory)
    {
        history->numValuesInXSum = 0;

        for (PullCoordinateHistory& pcrdHistory : history->pullCoordinateSums)
        {
            pcrdHistory.value    = 0;
            pcrdHistory.valueRef = 0;

            clear_dvec(pcrdHistory.dr01);
            clear_dvec(pcrdHistory.dr23);
            clear_dvec(pcrdHistory.dr45);
            clear_dvec(pcrdHistory.dynaX);
        }

        for (PullGroupHistory& pgrpHistory : history->pullGroupSums)
        {
            clear_dvec(pgrpHistory.x);
        }
    }
    if (resetFHistory)
    {
        history->numValuesInFSum = 0;
        for (PullCoordinateHistory& pcrdHistory : history->pullCoordinateSums)
        {
            pcrdHistory.scalarForce = 0;
        }
    }
}

static void pull_print_coord_dr_components(FILE* out, const ivec dim, const dvec dr, const int numValuesInSum)
{
    for (int m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", dr[m] / numValuesInSum);
        }
    }
}

template<typename T>
static void pull_print_coord_dr(FILE*                out,
                                const pull_params_t& pullParams,
                                const t_pull_coord&  coordParams,
                                const T&             pcrdData,
                                double               referenceValue,
                                const int            numValuesInSum)
{
    const double unit_factor = pull_conversion_factor_internal2userinput(coordParams);

    fprintf(out, "\t%g", pcrdData.value * unit_factor / numValuesInSum);

    if (pullParams.bPrintRefValue && coordParams.eType != PullingAlgorithm::External)
    {
        fprintf(out, "\t%g", referenceValue * unit_factor / numValuesInSum);
    }

    if (pullParams.bPrintComp)
    {
        pull_print_coord_dr_components(out, coordParams.dim, pcrdData.dr01, numValuesInSum);
        if (coordParams.ngroup >= 4)
        {
            pull_print_coord_dr_components(out, coordParams.dim, pcrdData.dr23, numValuesInSum);
        }
        if (coordParams.ngroup >= 6)
        {
            pull_print_coord_dr_components(out, coordParams.dim, pcrdData.dr45, numValuesInSum);
        }
    }
}

static void pull_print_x(FILE* out, pull_t* pull, double t)
{
    fprintf(out, "%.4f", t);

    for (size_t c = 0; c < pull->coord.size(); c++)
    {
        const pull_coord_work_t&     pcrd           = pull->coord[c];
        int                          numValuesInSum = 1;
        const PullCoordinateHistory* pcrdHistory    = nullptr;

        if (pull->bXOutAverage)
        {
            pcrdHistory = &pull->coordForceHistory->pullCoordinateSums[c];

            numValuesInSum = pull->coordForceHistory->numValuesInXSum;
            pull_print_coord_dr(
                    out, pull->params, pcrd.params_, *pcrdHistory, pcrdHistory->valueRef, numValuesInSum);
        }
        else
        {
            pull_print_coord_dr(
                    out, pull->params, pcrd.params_, pcrd.spatialData, pcrd.value_ref, numValuesInSum);
        }

        if (pull->params.bPrintCOM)
        {
            if (pcrd.params_.eGeom == PullGroupGeometry::Cylinder)
            {
                for (int m = 0; m < DIM; m++)
                {
                    if (pcrd.params_.dim[m])
                    {
                        /* This equates to if (pull->bXOutAverage) */
                        if (pcrdHistory)
                        {
                            fprintf(out, "\t%g", pcrdHistory->dynaX[m] / numValuesInSum);
                        }
                        else
                        {
                            fprintf(out, "\t%g", pcrd.dynamicGroup0->x[m]);
                        }
                    }
                }
            }
            else
            {
                for (int m = 0; m < DIM; m++)
                {
                    if (pcrd.params_.dim[m])
                    {
                        if (pull->bXOutAverage)
                        {
                            fprintf(out,
                                    "\t%g",
                                    pull->coordForceHistory->pullGroupSums[pcrd.params_.group[0]].x[m]
                                            / numValuesInSum);
                        }
                        else
                        {
                            fprintf(out, "\t%g", pull->group[pcrd.params_.group[0]].x[m]);
                        }
                    }
                }
            }
            for (int g = 1; g < pcrd.params_.ngroup; g++)
            {
                for (int m = 0; m < DIM; m++)
                {
                    if (pcrd.params_.dim[m])
                    {
                        if (pull->bXOutAverage)
                        {
                            fprintf(out,
                                    "\t%g",
                                    pull->coordForceHistory->pullGroupSums[pcrd.params_.group[g]].x[m]
                                            / numValuesInSum);
                        }
                        else
                        {
                            fprintf(out, "\t%g", pull->group[pcrd.params_.group[g]].x[m]);
                        }
                    }
                }
            }
        }
    }
    fprintf(out, "\n");

    if (pull->bXOutAverage)
    {
        pullResetHistory(pull->coordForceHistory, TRUE, FALSE);
    }
}

static void pull_print_f(FILE* out, const pull_t* pull, double t)
{
    fprintf(out, "%.4f", t);

    if (pull->bFOutAverage)
    {
        for (size_t c = 0; c < pull->coord.size(); c++)
        {
            fprintf(out,
                    "\t%g",
                    pull->coordForceHistory->pullCoordinateSums[c].scalarForce
                            / pull->coordForceHistory->numValuesInFSum);
        }
    }
    else
    {
        for (const pull_coord_work_t& coord : pull->coord)
        {
            fprintf(out, "\t%g", coord.scalarForce);
        }
    }
    fprintf(out, "\n");

    if (pull->bFOutAverage)
    {
        pullResetHistory(pull->coordForceHistory, FALSE, TRUE);
    }
}

void pull_print_output(struct pull_t* pull, int64_t step, double time)
{
    GMX_ASSERT(pull->numExternalPotentialsStillToBeAppliedThisStep == 0,
               "pull_print_output called before all external pull potentials have been applied");

    if (pull->params.nstxout != 0)
    {
        /* Do not update the average if the number of observations already equal (or are
         * higher than) what should be in each average output. This can happen when
         * appending to a file from a checkpoint, which would otherwise include the
         * last value twice.*/
        if (pull->bXOutAverage && !pull->coord.empty()
            && pull->coordForceHistory->numValuesInXSum < pull->params.nstxout)
        {
            addToPullxHistory(pull);
        }
        if (step % pull->params.nstxout == 0)
        {
            pull_print_x(pull->out_x, pull, time);
        }
    }

    if (pull->params.nstfout != 0)
    {
        /* Do not update the average if the number of observations already equal (or are
         * higher than) what should be in each average output. This can happen when
         * appending to a file from a checkpoint, which would otherwise include the
         * last value twice.*/
        if (pull->bFOutAverage && !pull->coord.empty()
            && pull->coordForceHistory->numValuesInFSum < pull->params.nstfout)
        {
            addToPullfHistory(pull);
        }
        if (step % pull->params.nstfout == 0)
        {
            pull_print_f(pull->out_f, pull, time);
        }
    }
}

static void set_legend_for_coord_components(const pull_coord_work_t*  pcrd,
                                            int                       coord_index,
                                            std::vector<std::string>* setname)
{
    /*  Loop over the distance vectors and print their components. Each vector is made up of two consecutive groups. */
    for (int g = 0; g < pcrd->params_.ngroup; g += 2)
    {
        /* Loop over the components */
        for (int m = 0; m < DIM; m++)
        {
            if (pcrd->params_.dim[m])
            {
                std::string legend;

                if (g == 0 && pcrd->params_.ngroup <= 2)
                {
                    /*  For the simplest case we print a simplified legend without group indices,
                       just the cooordinate index and which dimensional component it is. */
                    legend = gmx::formatString("%d d%c", coord_index + 1, 'X' + m);
                }
                else
                {
                    /* Otherwise, print also the group indices (relative to the coordinate, not the global ones). */
                    legend = gmx::formatString("%d g %d-%d d%c", coord_index + 1, g + 1, g + 2, 'X' + m);
                }

                setname->emplace_back(legend);
            }
        }
    }
}

static FILE* open_pull_out(const char*             fn,
                           struct pull_t*          pull,
                           const gmx_output_env_t* oenv,
                           gmx_bool                bCoord,
                           const bool              restartWithAppending)
{
    FILE* fp;
    int   m;

    std::vector<std::string> setname;
    if (restartWithAppending)
    {
        fp = gmx_fio_fopen(fn, "a+");
    }
    else
    {
        fp = gmx_fio_fopen(fn, "w+");
        if (bCoord)
        {
            auto buf = gmx::formatString("Position (nm%s)", pull->bAngle ? ", deg" : "");
            if (pull->bXOutAverage)
            {
                xvgr_header(fp, "Pull Average COM", "Time (ps)", buf, exvggtXNY, oenv);
            }
            else
            {
                xvgr_header(fp, "Pull COM", "Time (ps)", buf, exvggtXNY, oenv);
            }
        }
        else
        {
            auto buf = gmx::formatString("Force (kJ/mol/nm%s)", pull->bAngle ? ", kJ/mol/rad" : "");
            if (pull->bFOutAverage)
            {
                xvgr_header(fp, "Pull Average force", "Time (ps)", buf, exvggtXNY, oenv);
            }
            else
            {
                xvgr_header(fp, "Pull force", "Time (ps)", buf, exvggtXNY, oenv);
            }
        }

        /* With default mdp options only the actual coordinate value is printed (1),
         * but optionally the reference value (+ 1),
         * the group COMs for all the groups (+ ngroups_max*DIM)
         * and the components of the distance vectors can be printed (+ (ngroups_max/2)*DIM).
         */

        for (size_t c = 0; c < pull->coord.size(); c++)
        {
            if (bCoord)
            {
                /* The order of this legend should match the order of printing
                 * the data in print_pull_x above.
                 */

                /* The pull coord distance */
                setname.emplace_back(gmx::formatString("%zu", c + 1));
                if (pull->params.bPrintRefValue && pull->coord[c].params_.eType != PullingAlgorithm::External)
                {
                    setname.emplace_back(gmx::formatString("%zu ref", c + 1));
                }
                if (pull->params.bPrintComp)
                {
                    set_legend_for_coord_components(&pull->coord[c], c, &setname);
                }

                if (pull->params.bPrintCOM)
                {
                    for (int g = 0; g < pull->coord[c].params_.ngroup; g++)
                    {
                        /* Legend for reference group position */
                        for (m = 0; m < DIM; m++)
                        {
                            if (pull->coord[c].params_.dim[m])
                            {
                                setname.emplace_back(
                                        gmx::formatString("%zu g %d %c", c + 1, g + 1, 'X' + m));
                            }
                        }
                    }
                }
            }
            else
            {
                /* For the pull force we always only use one scalar */
                setname.emplace_back(gmx::formatString("%zu", c + 1));
            }
        }
        if (setname.size() > 1)
        {
            xvgrLegend(fp, setname, oenv);
        }
    }

    return fp;
}

void init_pull_output_files(pull_t*                     pull,
                            int                         nfile,
                            const t_filenm              fnm[],
                            const gmx_output_env_t*     oenv,
                            const gmx::StartingBehavior startingBehavior)
{
    /* Check for px and pf filename collision, if we are writing
       both files */
    std::string px_filename, pf_filename;
    std::string px_appended, pf_appended;
    try
    {
        px_filename = std::string(opt2fn("-px", nfile, fnm));
        pf_filename = std::string(opt2fn("-pf", nfile, fnm));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    bool restartWithAppending = startingBehavior == gmx::StartingBehavior::RestartWithAppending;
    if ((pull->params.nstxout != 0) && (pull->params.nstfout != 0) && (px_filename == pf_filename))
    {
        if (!opt2bSet("-px", nfile, fnm) && !opt2bSet("-pf", nfile, fnm))
        {
            /* We are writing both pull files but neither set directly. */
            try
            {
                px_appended = append_before_extension(px_filename, "_pullx");
                pf_appended = append_before_extension(pf_filename, "_pullf");
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            pull->out_x = open_pull_out(px_appended.c_str(), pull, oenv, TRUE, restartWithAppending);
            pull->out_f = open_pull_out(pf_appended.c_str(), pull, oenv, FALSE, restartWithAppending);
            return;
        }
        else
        {
            /* If at least one of -px and -pf is set but the filenames are identical: */
            gmx_fatal(FARGS, "Identical pull_x and pull_f output filenames %s", px_filename.c_str());
        }
    }
    if (pull->params.nstxout != 0)
    {
        pull->out_x = open_pull_out(opt2fn("-px", nfile, fnm), pull, oenv, TRUE, restartWithAppending);
    }
    if (pull->params.nstfout != 0)
    {
        pull->out_f = open_pull_out(opt2fn("-pf", nfile, fnm), pull, oenv, FALSE, restartWithAppending);
    }
}

void initPullHistory(pull_t* pull, ObservablesHistory* observablesHistory)
{
    GMX_RELEASE_ASSERT(pull, "Need a valid pull object");

    if (observablesHistory == nullptr)
    {
        pull->coordForceHistory = nullptr;
        return;
    }
    /* If pull->coordForceHistory is already set we are starting from a checkpoint. Do not reset it. */
    if (observablesHistory->pullHistory == nullptr)
    {
        observablesHistory->pullHistory          = std::make_unique<PullHistory>();
        pull->coordForceHistory                  = observablesHistory->pullHistory.get();
        pull->coordForceHistory->numValuesInXSum = 0;
        pull->coordForceHistory->numValuesInFSum = 0;
        pull->coordForceHistory->pullCoordinateSums.resize(pull->coord.size());
        pull->coordForceHistory->pullGroupSums.resize(pull->group.size());
    }
    else
    {
        pull->coordForceHistory = observablesHistory->pullHistory.get();
    }
}
