/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Defines the checkpoint data structure for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdtypes
 */

#include "gmxpre.h"

#include "checkpointdata.h"

#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/keyvaluetreeserializer.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

void ReadCheckpointDataHolder::deserialize(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(serializer->reading(),
                       "Tried to deserialize using a serializing ISerializer object.");

    checkpointTree_ = deserializeKeyValueTree(serializer);
}

void WriteCheckpointDataHolder::serialize(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(!serializer->reading(),
                       "Tried to serialize using a deserializing ISerializer object.");

    serializeKeyValueTree(outputTreeBuilder_.build(), serializer);

    // Tree builder should not be used after build() (see docstring)
    // Make new builder to leave object in valid state
    outputTreeBuilder_ = KeyValueTreeBuilder();
}

bool ReadCheckpointDataHolder::keyExists(const std::string& key) const
{
    return checkpointTree_.keyExists(key);
}

std::vector<std::string> ReadCheckpointDataHolder::keys() const
{
    std::vector<std::string> keys;
    for (const auto& property : checkpointTree_.properties())
    {
        keys.emplace_back(property.key());
    }
    return keys;
}

ReadCheckpointData ReadCheckpointDataHolder::checkpointData(const std::string& key) const
{
    return ReadCheckpointData(checkpointTree_[key].asObject());
}

void ReadCheckpointDataHolder::dump(FILE* out) const
{
    if (out != nullptr)
    {
        TextWriter textWriter(out);
        dumpKeyValueTree(&textWriter, checkpointTree_);
    }
}

WriteCheckpointData WriteCheckpointDataHolder::checkpointData(const std::string& key)
{
    hasCheckpointDataBeenRequested_ = true;
    return WriteCheckpointData(outputTreeBuilder_.rootObject().addObject(key));
}

bool WriteCheckpointDataHolder::empty() const
{
    return !hasCheckpointDataBeenRequested_;
}

} // namespace gmx
