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
#include <unordered_map>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

/*! \internal \brief
 * Neutron scattering factor parameters for an atom type.
 */
struct NeutronStructureFactor {

    //! atomic isotope
    std::string isotope;
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
 * X-ray scattering factor parameters for an atom type.
 */
struct XrayStructureFactor {
    //! atomic isotope
    std::string isotope;
    //! proton number
    int         proton;
    //! Parameters for the Cromer Mann fit
    CromerMan   CM;
};

//! Helper function to read in scattering data from file.
std::vector<NeutronStructureFactor> readNeutronScatteringFactors();

//! Helper function to read in scattering data from file.
std::vector<XrayStructureFactor> readXrayScatteringFactors();

/*! \internal \brief
 * Base class for computing x-ray and neutron scattering
 */
class ComputeScattering
{
    public:
        virtual ~ComputeScattering() = default;

        //! retrieves scattering length based on atom index
        virtual double getScatterLength(int i, double q) = 0;

        //! computes scattering intensity at zero scattering angle; for normalization
        double computeS0(Selection sel);

        //! computes the scattering intensity for a given scattering angle q
        double computeIntensity(t_pbc *pbc, double q, Selection sel);

        //! computes atomic form factors
        double getFormFactor(int i, int j, double q);
};

//! Helper to get the isotope of every atom
std::vector<std::string> getIsotopes(const t_atoms *atoms);

//! Hash function to allow use of pair in unordered_map
struct pairHash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>() (pair.first) ^ std::hash<T2>() (pair.second);
    }
};

/*! \internal \brief
 * Derrived class for computing x-ray scattering
 */
class XrayInfo : public ComputeScattering
{
    public:
        //! constructor
        XrayInfo(std::vector<std::string> isotopes, std::vector<double> qList);

        //! retrieves scattering length based on atom index
        double getScatterLength(int i, double q) override;

    private:
        //! refId of each atom in selection, should work for dynamic selections
        std::vector<std::string> isotopes_;

        //! scattering length of each atom in selection; same order as atomIndex_
        std::unordered_map<std::pair<std::string, double>, double, pairHash> scatterFactors_;

};

/*! \internal \brief
 * Derrived class for computing neutron scattering
 */
class NeutronInfo : public ComputeScattering
{
    public:
        //! constructor
        NeutronInfo(std::vector<std::string> isotopes);

        //! retrieves scattering length based on atom index
        double getScatterLength(int i, double q) override;

    private:
        //! refId of each atom in selection, should work for dynamic selections
        std::vector<std::string> isotopes_;

        //! scattering length of each atom in selection; same order as atomIndex_
        std::unordered_map<std::string, double> scatterFactors_;

};

} //namespace gmx
