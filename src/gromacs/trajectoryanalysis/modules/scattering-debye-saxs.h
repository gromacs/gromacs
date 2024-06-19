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
 * Declares class for SAXS Debye Scattering
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include <cstddef>

#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/topology/topology.h"

#include "isotope.h"
#include "scattering-debye.h"
#include "scatteringfactors.h"

namespace gmx
{

/*! \internal \brief
 * Hash function to allow use of pair in unordered_map
 */
struct pairHash
{
    template<class T1, class T2>

    //! there is no default operator so we define one
    std::size_t operator()(const std::pair<T1, T2>& pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};


/*! \internal \brief
 * Derrived class for computing neutron scattering
 */
class SaxsDebye : public ComputeDebyeScattering
{
public:
    //! constructor
    SaxsDebye(std::vector<Isotope> isotopes, const std::vector<double>& qList);
    //! retrieves scattering length based on atom index
    double getScatteringLength(int i, double q) override;

private:
    //! Vector containing enum of isotopes for each atom in selection
    std::vector<Isotope> isotopes_;

    //! scattering length of each atom in selection; same order as atomIndex_
    std::unordered_map<std::pair<int, double>, double, pairHash> scatterFactors_;
};

} // namespace gmx
