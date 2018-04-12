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

#include "gmxpre.h"

#include "scattering-forces.h"

#include <cmath>

#include <string>
#include <unordered_map>
#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

std::vector<NeutronStructureFactor> readNeutronScatteringFactors()
{
    std::vector<NeutronStructureFactor> neutronScatteringFactors;
    FILE *scatteringFile;
    scatteringFile = libopen("nsfactor.dat");
    char  line[1000];
    // all non-header lines
    while (get_a_line(scatteringFile, line, 1000))
    {
        int    proton;
        int    neutron;
        char   currentAtomType[8];
        double scatterLength;

        if (sscanf(line, "%s %d %d %lf", currentAtomType, &proton, &neutron, &scatterLength) == 4)
        {
            NeutronStructureFactor nsf = {currentAtomType, proton, neutron, scatterLength};
            neutronScatteringFactors.push_back(nsf);
        }
    }
    fclose(scatteringFile);
    return neutronScatteringFactors;
};

std::vector<XrayStructureFactor> readXrayScatteringFactors()
{
    std::vector<XrayStructureFactor> xrayScatteringFactors;
    FILE *scatteringFile;
    scatteringFile = libopen("sfactor.dat");
    char  line[1000];
    // all non-header lines
    while (get_a_line(scatteringFile, line, 1000))
    {
        int    proton;
        char   currentAtomType[8];
        double a1, a2, a3, a4, b1, b2, b3, b4, c;

        if (sscanf(line, "%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   currentAtomType, &proton, &a1, &a2, &a3, &a4,
                   &b1, &b2, &b3, &b4, &c) == 11)
        {
            std::vector<double> a   = {a1, a2, a3, a4};
            std::vector<double> b   = {b1, b2, b3, b4};
            CromerMan           CM  = {a, b, c};
            XrayStructureFactor xsf = {currentAtomType, proton, CM};
            xrayScatteringFactors.push_back(xsf);
        }
    }
    fclose(scatteringFile);
    return xrayScatteringFactors;
};

NeutronInfo::NeutronInfo(const t_atoms atoms)
{
    std::vector<NeutronStructureFactor>     neutronScatterFactors = readNeutronScatteringFactors();
    std::unordered_map<std::string, double> scatterFactors;
    for (auto scatter : neutronScatterFactors)
    {
        scatterFactors_[scatter.isotope] = scatter.scatterLength;
    }

    std::vector<std::string> isotopes;
    for (int i = 0; i < atoms.nr; ++i)
    {
        std::string isotope = atoms.atom[i].elem;
        isotopes_.push_back(isotope);
    }
};

XrayInfo::XrayInfo(const t_atoms atoms)
{
    std::vector<XrayStructureFactor>           xrayScatterFactors = readXrayScatteringFactors();
    std::unordered_map<std::string, CromerMan> scatterFactors;
    for (auto scatter : xrayScatterFactors)
    {
        scatterFactors_[scatter.isotope] = scatter.CM;
    }

    std::vector<std::string> isotopes;
    for (int i = 0; i < atoms.nr; ++i)
    {
        std::string isotope = atoms.atom[i].elem;
        isotopes_.push_back(isotope);
    }
}

double NeutronInfo::getScatterLength(int i, double gmx_unused q)
{
    std::string isotope    = isotopes_[i];
    double      scattering = scatterFactors_[isotope];
    return scattering;
};

double XrayInfo::getScatterLength(int i, double q)
{
    std::string isotope    = isotopes_[i];
    CromerMan   CM         = scatterFactors_[isotope];
    double      scattering = CM.c;
    // need factor of 10 to account for gromacs being in nm and
    // scatter factors being in Å
    double q4pi = 10 * q / (4 * M_PI);
    for (int j = 0; j < 4; j++)
    {
        scattering += CM.a[j] * exp(-CM.b[j] * q4pi * q4pi);
    }
    return scattering;
};

double ComputeScattering::computeS0(Selection sel)
{
    //significant adding errors if Debeye is not a double
    double    DebeyeS0 = 0;
    const int posCount = sel.posCount();
    for (int i = 0; i < posCount; ++i)
    {
        const int index_i = sel.position(i).atomIndices()[0];
        for (int j = i + 1; j < posCount; ++j)
        {
            const int index_j = sel.position(j).atomIndices()[0];
            DebeyeS0 += getFormFactor(index_i, index_j, 0);
        }
    }
    return DebeyeS0;
};

double ComputeScattering::computeIntensity(t_pbc *pbc, double q, Selection sel)
{
    //significant adding errors if Debeye is not a double
    double    Debeye   = 0;
    const int posCount = sel.posCount();
    for (int i = 0; i < posCount; ++i)
    {
        const SelectionPosition &pi      = sel.position(i);
        const int                index_i = sel.position(i).atomIndices()[0];
        for (int j = i + 1; j < posCount; ++j)
        {
            const SelectionPosition &pj         = sel.position(j);
            const int                index_j    = sel.position(j).atomIndices()[0];
            double                   formFactor = getFormFactor(index_i, index_j, q);
            rvec                     dx;
            if (pbc != nullptr)
            {
                pbc_dx(pbc, pi.x(), pj.x(), dx);
            }
            else
            {
                rvec_sub(pi.x(), pj.x(), dx);
            }
            double qDist    = q * std::sqrt(norm2(dx));
            double sinQDist = sin(qDist);
            Debeye += formFactor * (sinQDist / qDist);
        }
    }
    return Debeye;
};

double ComputeScattering::getFormFactor(int i, int j, double q)
{
    auto scatterLength_i = getScatterLength(i, q);
    auto scatterLength_j = getScatterLength(j, q);
    // multiply by 2 since we only store half the form factors (ij == ji)
    return 2 * scatterLength_i * scatterLength_j;
}

}; //namespace gmx
