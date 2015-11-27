/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief Implementions of related classes for tests that want to
 * inspect energies produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "energyreader.h"

#include <cstdio>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

EnergyFileReader::EnergyFileReader(const std::string &filename)
    : filename_(filename),
      energyFile_(nullptr),
      frame_(),
      frameContents_(&frame_)
{
    init_enxframe(&frame_);
}

EnergyFrameReader
EnergyFileReader::openToReadFields(const std::vector<std::string> &namesOfRequiredEnergyFields)
{
    if (energyFile_)
    {
        GMX_THROW(APIError("Opening energy file " + filename_ + " twice is not supported"));
    }
    energyFile_.reset(open_enx(filename_.c_str(), "r"));
    if (!energyFile_)
    {
        GMX_THROW(FileIOError("Could not open energy file " + filename_ + " for reading"));
    }

    /* Read in the names of energy fields used in this file. The
     * resulting data structure would leak if an exception was thrown,
     * so transfer the contents that we actually need to a map we can
     * keep.
     *
     * TODO Technically, the insertions into the map could throw
     * std::bad_alloc and we could leak memory allocated by
     * do_enxnms(), but there's nothing we can do about this right
     * now. */
    std::map<std::string, int> indicesOfEnergyFields;
    {
        int          numEnergyTerms;
        gmx_enxnm_t *energyNames = nullptr;
        do_enxnms(energyFile_.get(), &numEnergyTerms, &energyNames);
        for (int i = 0; i != numEnergyTerms; ++i)
        {
            const char *name           = energyNames[i].name;
            auto        requiredEnergy = std::find_if(std::begin(namesOfRequiredEnergyFields),
                                                      std::end(namesOfRequiredEnergyFields),
                                                      [name](const std::string &n){
                                                          return 0 == n.compare(name);
                                                      });
            if (requiredEnergy != namesOfRequiredEnergyFields.end())
            {
                indicesOfEnergyFields[name] = i;
            }
        }
        // Clean up old data structures
        free_enxnms(numEnergyTerms, energyNames);
    }

    // Throw if we failed to find the fields we need
    if (indicesOfEnergyFields.size() != namesOfRequiredEnergyFields.size())
    {
        std::string requiredEnergiesNotFound = "Did not find the following required energies in mdrun output:\n";
        for (auto &name : namesOfRequiredEnergyFields)
        {
            auto possibleIndex = indicesOfEnergyFields.find(name);
            if (possibleIndex == indicesOfEnergyFields.end())
            {
                requiredEnergiesNotFound += name + "\n";
            }
        }
        GMX_THROW(APIError(requiredEnergiesNotFound));
    }

    return EnergyFrameReader(indicesOfEnergyFields, energyFile_.get(), &frame_);
}

EnergyFrameReader::EnergyFrameReader(std::map<std::string, int> indicesOfEnergyFields,
                                     ener_file *energyFile,
                                     t_enxframe *frame)
    : indicesOfEnergyFields_(indicesOfEnergyFields),
      energyFile_(energyFile),
      frame_(frame)
{
}

EnergyFrameInfo
EnergyFrameReader::readNextFrame()
{
    EnergyFrameInfo frameInfo;

    bool            foundFrame = do_enx(energyFile_, frame_);
    // Hack to help the output here cope with how the .edr-reading code writes output
    std::fputc('\n', stderr);
    if (!foundFrame)
    {
        return frameInfo;
    }
    frameInfo.time_    = frame_->t;
    frameInfo.step_    = frame_->step;
    frameInfo.isValid_ = true;

    // Copy the energy value that was read into the EnergyFrameInfo to return it
    for (auto &index : indicesOfEnergyFields_)
    {
        if (index.second >= frame_->nre)
        {
            GMX_THROW(InternalError(formatString("Index %d for energy %s not present in energy frame with %d energies",
                                                 index.second, index.first.c_str(), frame_->nre)));
        }
        frameInfo.values_[index.first] = frame_->ener[index.second].e;
    }
    return frameInfo;
}

EnergyFrameInfo::operator bool() const
{
    return isValid_;
}

std::string EnergyFrameInfo::getFrameName() const
{
    if (!isValid_)
    {
        GMX_THROW(APIError("Cannot get energy-frame info unless a valid frame has been read"));
    }
    return formatString("Time %f Step %" GMX_PRId64, time_, step_);
}

real EnergyFrameInfo::getValue(const std::string &name) const
{
    if (!isValid_)
    {
        GMX_THROW(APIError("Cannot get energy-frame information unless a valid frame has been read"));
    }
    auto valueIterator = values_.find(name);
    if (valueIterator != values_.end())
    {
        return valueIterator->second;
    }
    GMX_THROW(APIError("Cannot get energy value " + name + " unless previously registered when constructing EnergyFrameReader"));
}

} // namespace
} // namespace
