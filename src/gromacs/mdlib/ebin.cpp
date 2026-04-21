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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "ebin.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <filesystem>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"

int t_ebin::getSpace(gmx::ArrayRef<const std::string> enm, const char* unit)
{
    const int index = nener_;
    nener_ += enm.ssize();
    accumulation_.energies_.resize(nener_);
    simulationAccumulation_.energies_.resize(nener_);
    enm_.resize(nener_);
    for (int i = index; i < nener_; i++)
    {
        enm_[i].name = enm[i - index];
        if (unit != nullptr)
        {
            enm_[i].unit = unit;
        }
        else
        {
            /* Determine the unit from the longname.
             * These units should have been defined in ifunc.c
             * But even better would be if all interactions functions
             * return energies and all non-interaction function
             * entries would be removed from the ifunc array.
             */
            const char* u = unit_energy;
            for (const auto f : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                if (std::strcmp(enm_[i].name.c_str(), interaction_function[f].longname) == 0)
                {
                    /* Only the terms in this list are not energies */
                    switch (f)
                    {
                        case InteractionFunction::DistanceRestraintViolations:
                            u = unit_length;
                            break;
                        case InteractionFunction::OrientationRestraintDeviations: u = "obs"; break;
                        case InteractionFunction::Temperature: u = unit_temp_K; break;
                        case InteractionFunction::PressureDispersionCorrection:
                        case InteractionFunction::Pressure: u = unit_pres_bar; break;
                        default: break;
                    }
                }
            }
            enm_[i].unit = u;
        }
    }

    return index;
}

void t_ebin::addValues(const int entryIndex, const gmx::ArrayRef<const real> ener, const bool accumulate)
{
    if ((entryIndex + gmx::ssize(ener) > nener_) || (entryIndex < 0))
    {
        gmx_fatal(FARGS,
                  "%s-%d: Energies out of range: entryIndex=%d nener=%d maxener=%d",
                  __FILE__,
                  __LINE__,
                  entryIndex,
                  static_cast<int>(gmx::ssize(ener)),
                  nener_);
    }

    for (int i = 0; i < gmx::ssize(ener); i++)
    {
        accumulation_.energies_[entryIndex + i].e = ener[i];
    }

    if (accumulate)
    {
        const int m = accumulation_.sumCount_;

        if (m == 0)
        {
            for (int i = 0; i < gmx::ssize(ener); i++)
            {
                accumulation_.energies_[entryIndex + i].sumSqDev = 0;
                accumulation_.energies_[entryIndex + i].esum     = ener[i];
                simulationAccumulation_.energies_[entryIndex + i].esum += ener[i];
            }
        }
        else
        {
            const double invmm = (1.0 / m) / (m + 1.0);

            for (int i = 0; i < gmx::ssize(ener); i++)
            {
                /* Value for this component */
                const double e = ener[i];

                /* first update sigma, then sum */
                const double diff = accumulation_.energies_[entryIndex + i].esum - m * e;

                accumulation_.energies_[entryIndex + i].sumSqDev += diff * diff * invmm;
                accumulation_.energies_[entryIndex + i].esum += e;
                simulationAccumulation_.energies_[entryIndex + i].esum += e;
            }
        }
    }
}

void t_ebin::addValue(int entryIndex, const real ener, bool accumulate)
{
    addValues(entryIndex, gmx::constArrayRefFromArray(&ener, 1), accumulate);
}

// TODO It would be faster if this function was templated on both bSum
// and whether eb->nsum was zero, to lift the branches out of the loop
// over all possible energy terms, but that is true for a lot of the
// ebin and mdebin functionality, so we should do it all or nothing.
void t_ebin::addValuesIndexed(const int                 entryIndex,
                              gmx::ArrayRef<const bool> shouldUse,
                              gmx::ArrayRef<const real> ener,
                              const bool                accumulate)
{
    GMX_ASSERT(shouldUse.size() == ener.size(), "View sizes must match");
    GMX_ASSERT(entryIndex + std::count(shouldUse.begin(), shouldUse.end(), true) <= nener_,
               gmx::formatString("Energies out of range: entryIndex=%d nener=%td maxener=%d",
                                 entryIndex,
                                 std::count(shouldUse.begin(), shouldUse.end(), true),
                                 nener_)
                       .c_str());
    GMX_ASSERT(entryIndex >= 0, "Must have non-negative entry");

    const int    m              = accumulation_.sumCount_;
    const double invmm          = (m == 0) ? 0 : (1.0 / m) / (m + 1.0);
    t_energy*    energyEntry    = &accumulation_.energies_[entryIndex];
    t_energy*    simEnergyEntry = &simulationAccumulation_.energies_[entryIndex];
    auto         shouldUseIter  = shouldUse.begin();
    for (const auto& theEnergy : ener)
    {
        if (*shouldUseIter)
        {
            energyEntry->e = theEnergy;
            if (accumulate)
            {
                if (m == 0)
                {
                    energyEntry->sumSqDev = 0;
                    energyEntry->esum     = theEnergy;
                    simEnergyEntry->esum += theEnergy;
                }
                else
                {
                    /* first update sigma, then sum */
                    double diff = energyEntry->esum - m * theEnergy;
                    energyEntry->sumSqDev += diff * diff * invmm;
                    energyEntry->esum += theEnergy;
                    simEnergyEntry->esum += theEnergy;
                }
                ++simEnergyEntry;
            }
            ++energyEntry;
        }
        ++shouldUseIter;
    }
}

void t_ebin::incrementCount(const bool increaseSums)
{
    accumulation_.numSteps_ += 1;
    simulationAccumulation_.numSteps_ += 1;

    if (increaseSums)
    {
        accumulation_.sumCount_ += 1;
        simulationAccumulation_.sumCount_ += 1;
    }
}

void t_ebin::resetSums()
{
    accumulation_.numSteps_ = 0;
    accumulation_.sumCount_ = 0;
    /* The actual sums are cleared when the next frame is stored */
}

void t_ebin::restoreFromEnergyHistory(const energyhistory_t& enerhist)
{
    GMX_RELEASE_ASSERT((enerhist.nsum == 0 || nener_ == gmx::ssize(enerhist.ener_sum))
                               && (enerhist.nsum_sim == 0 || nener_ == gmx::ssize(enerhist.ener_sum_sim)),
                       "Items counts in history should match the simulation");

    accumulation_.numSteps_           = enerhist.nsteps;
    accumulation_.sumCount_           = enerhist.nsum;
    simulationAccumulation_.numSteps_ = enerhist.nsteps_sim;
    simulationAccumulation_.sumCount_ = enerhist.nsum_sim;

    for (int i = 0; i < nener_; i++)
    {
        accumulation_.energies_[i].sumSqDev = (enerhist.nsum > 0 ? enerhist.ener_sumSqDev[i] : 0);
        accumulation_.energies_[i].esum     = (enerhist.nsum > 0 ? enerhist.ener_sum[i] : 0);
        simulationAccumulation_.energies_[i].esum =
                (enerhist.nsum_sim > 0 ? enerhist.ener_sum_sim[i] : 0);
    }
}

void pr_ebin(FILE* fp, const t_ebin& eb, int entryIndex, int nener, int nperline, int prmode, bool bPrHead)
{
    int  i, j, i0;
    int  rc;
    char buf[30];

    rc = 0;

    if (entryIndex < 0 || entryIndex > eb.numTerms())
    {
        gmx_fatal(FARGS, "Invalid entryIndex in pr_ebin: %d", entryIndex);
    }
    int start = entryIndex;
    if (nener > eb.numTerms())
    {
        gmx_fatal(FARGS, "Invalid nener in pr_ebin: %d", nener);
    }
    int end = eb.numTerms();
    if (nener != -1)
    {
        end = entryIndex + nener;
    }
    for (i = start; (i < end) && rc >= 0;)
    {
        if (bPrHead)
        {
            i0 = i;
            for (j = 0; (j < nperline) && (i < end) && rc >= 0; j++, i++)
            {
                if (std::strncmp(eb.names()[i].name.c_str(), "Pres", 4) == 0)
                {
                    /* Print the pressure unit to avoid confusion */
                    sprintf(buf, "%s (%s)", eb.names()[i].name.c_str(), unit_pres_bar);
                    rc = fprintf(fp, "%15s", buf);
                }
                else
                {
                    rc = fprintf(fp, "%15s", eb.names()[i].name.c_str());
                }
            }

            if (rc >= 0)
            {
                rc = fprintf(fp, "\n");
            }

            i = i0;
        }
        for (j = 0; (j < nperline) && (i < end) && rc >= 0; j++, i++)
        {
            switch (prmode)
            {
                case eprNORMAL:
                    rc = fprintf(fp, "   %12.5e", eb.accumulation().energies()[i].e);
                    break;
                case eprAVER:
                    if (eb.simulationAccumulation().sumCount() > 0)
                    {
                        rc = fprintf(fp, "   %12.5e", eb.simulationAccumulation().average(i));
                    }
                    else
                    {
                        rc = fprintf(fp, "    %-12s", "N/A");
                    }
                    break;
                default: gmx_fatal(FARGS, "Invalid print mode %d in pr_ebin", prmode);
            }
        }
        if (rc >= 0)
        {
            rc = fprintf(fp, "\n");
        }
    }
    if (rc < 0)
    {
        gmx_fatal(FARGS, "Cannot write to logfile; maybe you are out of disk space?");
    }
}
