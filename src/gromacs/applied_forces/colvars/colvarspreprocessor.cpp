/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * \brief
 * Implements the Colvars GROMACS proxy class during pre-processing.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "colvarspreprocessor.h"

#include <list>
#include <map>
#include <sstream>
#include <string>

#include "gromacs/applied_forces/colvars/colvarproxygromacs.h"
#include "gromacs/pbcutil/pbc.h"

enum class PbcType : int;

namespace gmx
{
class MDLogger;

ColvarsPreProcessor::ColvarsPreProcessor(const std::string&   colvarsConfigString,
                                         t_atoms              atoms,
                                         PbcType              pbcType,
                                         const MDLogger*      logger,
                                         real                 ensembleTemperature,
                                         int                  seed,
                                         const matrix         box,
                                         ArrayRef<const RVec> x) :
    ColvarProxyGromacs(colvarsConfigString,
                       atoms,
                       pbcType,
                       logger,
                       true,
                       std::map<std::string, std::string>(),
                       ensembleTemperature,
                       seed),
    x_(x)
{

    // Initialize t_pbc struct
    set_pbc(&gmxPbc_, pbcType, box);

    cvm::log(cvm::line_marker);
    cvm::log("End colvars Initialization.\n\n");
}

std::vector<RVec> ColvarsPreProcessor::getColvarsCoords()
{

    std::vector<RVec> colvarsCoords;
    colvarsCoords.reserve(atoms_ids.size());
    for (const auto& atom_id : atoms_ids)
    {
        colvarsCoords.push_back(x_[atom_id]);
    }
    return colvarsCoords;
}

bool ColvarsPreProcessor::inputStreamsToKVT(KeyValueTreeObjectBuilder treeBuilder, const std::string& tag)
{

    // Save full copy of the content of the input streams (aka input files) into the KVT.
    for (const auto& inputName : list_input_stream_names())
    {
        std::istream&      stream = input_stream(inputName);
        std::ostringstream os;
        os << stream.rdbuf();
        std::string key = tag;
        key += "-";
        key += inputName;
        treeBuilder.addValue<std::string>(key, os.str());
    }
    return true;
}


} // namespace gmx
