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

#include <string>
#include <vector>
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

namespace gmx
{

/*! \internal \brief
 * Base class for storing scattering factor parameters for an atom type.
 */
typedef struct neutronStructureFactor {
    int          p;              /* proton number */
    int          n;              /* neutron number */
    real         scatterLength;  /* neutron scattering length*/
    std::string  element;        /* atomic element */
} neutronStructureFactor;

/*! \internal \brief
 * Helper function to read in scattering data from file.
 * TODO: This should really be a member of the ComputeScattering class but I 
 *       could not figure out how to make it work with the factory function
 */
void readScatteringFactors(std::vector<neutronStructureFactor> &neutronScatteringFactors_);

class ComputeScattering
{

private:
    //! refId of each atom in selection, should work for dynamic selections
    std::vector<int> atomIndex_;

    //! scattering length of each atom in selection; same order as atomIndex_
    std::vector<real> scatteringLength_;

    //! force on each atom in scattering group; same order as atomIndex_
    std::vector<real> scatterForces_;

    //! force constant used for calculationg forces
    real forceConstant_ = 1;

public:
    ComputeScattering(std::vector<int> &atomIndex_, std::vector<real> scatteringLength_) :
        atomIndex_(atomIndex_), scatteringLength_(scatteringLength_)  {};
    friend ComputeScattering makeScattering(const Selection &sel, const TopologyInformation &top);

    //! computes scattering intensity at zero scattering angle; for normalization
    double computeS0();

    //! retrieves scattering length based on selection refId
    real getScatterLength(int i);
    
    //! computes the scattering intensity for a given scattering angle q
    //! TODO: make the ComputeScattering class have sel as a member so it
    //! does not need to be passed as an input
    double computeIntensity(t_pbc *pbc, double q, Selection sel);

    void computeForces(t_pbc *pbc, double q, Selection sel, double scatterDiff);

};

//! factory function that returns an initialized ComputeScattering object
ComputeScattering makeScattering(const Selection &sel, const TopologyInformation &top);

} //namespace gmx
