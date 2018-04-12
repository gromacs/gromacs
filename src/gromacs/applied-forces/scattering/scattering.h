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

#include <map>
#include <string>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

namespace gmx
{

/*! \internal \brief
 * Struct for storing neutron scattering factor parameters for an atom type.
 */
struct neutronStructureFactor {

    //! atomic element
    std::string element;
    //! proton number
    int         proton;
    //! neutron number
    int         neutron;
    //! neutron scattering length
    double      scatterLength;
};

/*! \internal \brief
 * Cromer-Man scattering factor parameters for an atom type.
 */
struct CromerMan {
    //! parameter a
    std::vector<double> a;
    //! parameter a
    std::vector<double> b;
    //! parameter a
    double              c;
};

/*! \internal \brief
 * Struct for storing x-ray scattering factor parameters for an atom type.
 */
struct xrayStructureFactor {
    //! atomic element
    std::string element;
    //! proton number
    int         proton;
    //! Parameters for the Cromer Mann fit
    CromerMan   CM;
};

/*! \internal \brief
 * Helper function to read in scattering data from file.
 */
void readNeutronScatteringFactors(std::vector<neutronStructureFactor> &neutronScatteringFactors_);

/*! \internal \brief
 * Helper function to read in scattering data from file.
 */
void readXrayScatteringFactors(std::vector<xrayStructureFactor> &xrayScatteringFactors_);

class ComputeScattering
{
    public:
        virtual ~ComputeScattering() = default;

        //ComputeScattering() {std::move(*this);}

        //! computes scattering intensity at zero scattering angle; for normalization
        virtual double computeS0(Selection sel) = 0;

        //! computes the scattering intensity for a given scattering angle q
        virtual double computeIntensity(t_pbc *pbc, double q, Selection sel) = 0;
};

class XrayInfo : public ComputeScattering
{
    private:
        //! refId of each atom in selection, should work for dynamic selections
        std::vector<std::string> elements_;

        //! scattering length of each atom in selection; same order as atomIndex_
        std::map<std::string, CromerMan> scatterFactors_;

        //! retrieves scattering length based on atom index
        double getScatterLength(int i, double q);

    public:
        XrayInfo(std::vector<std::string> elements,
                 std::map<std::string, CromerMan> scatterFactors) :
            elements_(elements),
            scatterFactors_(scatterFactors)
        {
        }

        //! computes atomic form factors
        double getFormFactor(int i, int j, double q);

        //! computes scattering intensity at zero scattering angle; for normalization
        double computeS0(Selection sel);

        //! computes the scattering intensity for a given scattering angle q
        double computeIntensity(t_pbc *pbc, double q, Selection sel);
};

class NeutronInfo : public ComputeScattering
{
    private:
        //! refId of each atom in selection, should work for dynamic selections
        std::vector<std::string> elements_;

        //! scattering length of each atom in selection; same order as atomIndex_
        std::map<std::string, double> scatterFactors_;

        //! retrieves scattering length based on atom index
        double getScatterLength(int i);

    public:
        NeutronInfo(std::vector<std::string> elements,
                    std::map<std::string, double> scatterFactors) :
            elements_(elements),
            scatterFactors_(scatterFactors)
        {
        }

        //! computes atomic form factors
        double getFormFactor(int i, int j, double q);

        //! computes scattering intensity at zero scattering angle; for normalization
        double computeS0(Selection sel);

        //! computes the scattering intensity for a given scattering angle q
        double computeIntensity(t_pbc *pbc, double q, Selection sel);
};

//! factory function that returns an initialized NeutronInfo class
NeutronInfo makeNeutronInfo(const TopologyInformation &top);

//! factory function that returns an initialized XrayInfo class
XrayInfo makeXrayInfo(const TopologyInformation &top);

} //namespace gmx
