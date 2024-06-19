/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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

/*! \internal \file
 * \brief Implementions of related classes for tests that want to
 * inspect energies produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "energyreader.h"

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/fileio/enxio.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

struct ener_file;

namespace gmx
{
namespace test
{

EnergyFrameReaderPtr openEnergyFileToReadTerms(const std::string&              filename,
                                               const std::vector<std::string>& namesOfRequiredEnergyTerms)
{
    ener_file_ptr energyFile(open_enx(filename.c_str(), "r"));

    if (!energyFile)
    {
        GMX_THROW(FileIOError("Could not open energy file " + filename + " for reading"));
    }

    /* Read in the names of energy terms used in this file. The
     * resulting data structure would leak if an exception was thrown,
     * so transfer the contents that we actually need to a map we can
     * keep.
     *
     * TODO Technically, the insertions into the map could throw
     * std::bad_alloc and we could leak memory allocated by
     * do_enxnms(), but there's nothing we can do about this right
     * now. */
    std::map<std::string, int> indicesOfEnergyTerms;
    {
        int          numEnergyTerms;
        gmx_enxnm_t* energyNames = nullptr;
        do_enxnms(energyFile.get(), &numEnergyTerms, &energyNames);
        for (int i = 0; i != numEnergyTerms; ++i)
        {
            const char* name           = energyNames[i].name;
            auto        requiredEnergy = std::find_if(std::begin(namesOfRequiredEnergyTerms),
                                               std::end(namesOfRequiredEnergyTerms),
                                               [name](const std::string& n) { return name == n; });
            if (requiredEnergy != namesOfRequiredEnergyTerms.end())
            {
                indicesOfEnergyTerms[name] = i;
            }
        }
        // Clean up old data structures
        free_enxnms(numEnergyTerms, energyNames);
    }

    // Throw if we failed to find the terms we need
    if (indicesOfEnergyTerms.size() != namesOfRequiredEnergyTerms.size())
    {
        std::string requiredEnergiesNotFound =
                "Did not find the following required energies in mdrun output:\n";
        for (const auto& name : namesOfRequiredEnergyTerms)
        {
            auto possibleIndex = indicesOfEnergyTerms.find(name);
            if (possibleIndex == indicesOfEnergyTerms.end())
            {
                requiredEnergiesNotFound += name + "\n";
            }
        }
        GMX_THROW(APIError(requiredEnergiesNotFound));
    }

    return EnergyFrameReaderPtr(
            std::make_unique<EnergyFrameReader>(indicesOfEnergyTerms, energyFile.release()));
}

//! Helper function to obtain resources
static t_enxframe* make_enxframe()
{
    t_enxframe* frame;

    snew(frame, 1);
    init_enxframe(frame);

    return frame;
}

//! Helper function to clean up resources
void done_enxframe(t_enxframe* fr)
{
    // Free the contents, then the pointer itself
    free_enxframe(fr);
    sfree(fr);
}

// === EnergyFrameReader ===

EnergyFrameReader::EnergyFrameReader(const std::map<std::string, int>& indicesOfEnergyTerms,
                                     ener_file*                        energyFile) :
    indicesOfEnergyTerms_(indicesOfEnergyTerms),
    energyFileGuard_(energyFile),
    enxframeGuard_(make_enxframe()),
    haveProbedForNextFrame_(false),
    nextFrameExists_(false)
{
}

bool EnergyFrameReader::readNextFrame()
{
    if (haveProbedForNextFrame_)
    {
        if (nextFrameExists_)
        {
            GMX_THROW(
                    APIError("This frame has already been probed for, it should be used before "
                             "probing again."));
        }
        else
        {
            GMX_THROW(
                    APIError("This frame has already been probed for, it doesn't exist, so there "
                             "should not be subsequent attempts to probe for it."));
        }
    }
    haveProbedForNextFrame_ = true;
    // If there's a next frame, read it into enxframe_, and report the result.
    return nextFrameExists_ = do_enx(energyFileGuard_.get(), enxframeGuard_.get());
}

EnergyFrame EnergyFrameReader::frame()
{
    if (!haveProbedForNextFrame_)
    {
        readNextFrame();
    }
    if (!nextFrameExists_)
    {
        GMX_THROW(
                APIError("There is no next frame, so there should have been no attempt to use the "
                         "data, e.g. by reacting to a call to readNextFrame()."));
    }

    // Prepare for reading future frames
    haveProbedForNextFrame_ = false;
    nextFrameExists_        = false;

    // The probe filled enxframe_ with new data, so now we use that data to fill energyFrame
    return EnergyFrame(*enxframeGuard_.get(), indicesOfEnergyTerms_);
}

} // namespace test
} // namespace gmx
