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
#include <vector>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

std::vector<NeutronStructureFactor> readNeutronScatteringFactors()
{
    std::vector<NeutronStructureFactor> neutronScatteringFactors;
    gmx::FilePtr fp = gmx::openLibraryFile("nsfactor.dat");
    char         line[1000];
    // all non-header lines
    while (get_a_line(fp.get(), line, 1000))
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
    return neutronScatteringFactors;
};

std::vector<XrayStructureFactor> readXrayScatteringFactors()
{
    std::vector<XrayStructureFactor> xrayScatteringFactors;
    gmx::FilePtr                     fp = gmx::openLibraryFile("sfactor.dat");
    char  line[1000];
    // all non-header lines
    while (get_a_line(fp.get(), line, 1000))
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
    return xrayScatteringFactors;
};

std::vector<std::string> getIsotopes(const t_atoms *atoms)
{
    std::vector<std::string> isotopes;
    for (int i = 0; i < atoms->nr; ++i)
    {
        std::string isotope = atoms->atom[i].elem;
        isotopes.push_back(isotope);
    }
    return isotopes;
};

NeutronInfo::NeutronInfo(std::vector<std::string> isotopes)
{
    std::vector<NeutronStructureFactor>     neutronScatterFactors = readNeutronScatteringFactors();
    for (auto scatter : neutronScatterFactors)
    {
        scatterFactors_[scatter.isotope] = scatter.scatterLength;
    }
    isotopes_ = std::move(isotopes);
};

XrayInfo::XrayInfo(std::vector<std::string> isotopes, std::vector<double> qList)
{
    std::vector<XrayStructureFactor>           xrayScatterFactors = readXrayScatteringFactors();
    for (auto q : qList)
    {
        for (auto scatter : xrayScatterFactors)
        {

            double      scattering = scatter.CM.c;
            // need factor of 10 to convert from nm to Å
            double      q4pi = NM2A * q / (4 * M_PI);
            for (int j = 0; j < 4; j++)
            {
                scattering += scatter.CM.a[j] * exp(-scatter.CM.b[j] * q4pi * q4pi);
            }
            auto pair = std::make_pair(scatter.isotope, q);
            scatterFactors_[pair] = scattering;
        }
    }
    isotopes_ = std::move(isotopes);
}

double NeutronInfo::getScatterLength(int i, double gmx_unused q)
{
    std::string isotope    = isotopes_[i];
    double      scattering = scatterFactors_[isotope];
    return scattering;
};

double XrayInfo::getScatterLength(int i, double q)
{
    auto   pair       = std::make_pair(isotopes_[i], q);
    double scattering = scatterFactors_[pair];
    return scattering;
};

double ComputeScattering::computeS0MC(Selection sel, float coverage, int seed)
{
    //significant adding errors if Debeye is not a double
    double                      DebeyeS0 = 0;
    const int                   posCount = sel.posCount();
    DefaultRandomEngine         rng(seed);
    UniformIntDistribution<int> distribution(0, posCount-1);
    int                         scaledPosCount = int(coverage * posCount);
    if (scaledPosCount <= 0 || scaledPosCount > posCount)
    {
        gmx_fatal(FARGS, "Failed!\n");
    }
    for (int i = 0; i < scaledPosCount; ++i)
    {
        int       rand_i  = distribution(rng);
        const int index_i = sel.position(rand_i).atomIndices()[0];
        for (int j = 0; j < scaledPosCount; ++j)
        {
            int       rand_j = distribution(rng);
            if (rand_i == rand_j)
            {
                continue;
            }
            const int index_j = sel.position(rand_j).atomIndices()[0];
            DebeyeS0 += getFormFactor(index_i, index_j, 0);
        }
    }
    return DebeyeS0;
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
            // multiply by 2 since we only store half the form factors (ij == ji)
            DebeyeS0 += 2 * getFormFactor(index_i, index_j, 0);
        }
    }
    return DebeyeS0;
};

double ComputeScattering::computeIntensityMC(t_pbc *pbc, double q, Selection sel, float coverage, int seed)
{
    //significant adding errors if Debeye is not a double
    double                      Debeye   = 0;
    const int                   posCount = sel.posCount();
    DefaultRandomEngine         rng(seed);
    UniformIntDistribution<int> distribution(0, posCount-1);
    int                         scaledPosCount = int(coverage * posCount);
    if (scaledPosCount <= 0 || scaledPosCount > posCount)
    {
        gmx_fatal(FARGS, "You wound up with more points to compute than there were atoms. Try a lower coverage.\n");
    }
    for (int i = 0; i < scaledPosCount; ++i)
    {
        int                      rand_i  = distribution(rng);
        const SelectionPosition &pi      = sel.position(rand_i);
        const int                index_i = sel.position(rand_i).atomIndices()[0];
        for (int j = 0; j < scaledPosCount; ++j)
        {
            int                      rand_j = distribution(rng);
            if (rand_i == rand_j)
            {
                continue;
            }
            const SelectionPosition &pj         = sel.position(rand_j);
            const int                index_j    = sel.position(rand_j).atomIndices()[0];
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
            // multiply by 2 since we only store half the form factors (ij == ji)
            double                   formFactor = 2 * getFormFactor(index_i, index_j, q);
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
    return scatterLength_i * scatterLength_j;
}

}; //namespace gmx
