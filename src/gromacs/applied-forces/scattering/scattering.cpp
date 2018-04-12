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


#include "scattering.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/topology/topology.h"

#include <map>
#include <string>
#include <vector>

namespace gmx
{

void readScatteringFactors(std::vector<neutronStructureFactor> &neutronScatteringFactors_)
{
    FILE *scatteringFile;
    scatteringFile = libopen("nsfactor.dat");
    std::map < std::string, neutronStructureFactor> neutronStructureFactorLookup;
    char line[1000];
    // all non-header lines
    while (get_a_line(scatteringFile, line, 1000))
    {

        char  currentAtomType[8];
        int   p;
        int   n;
        float  scatterLength;

        if (sscanf(line, "%s %d %d %f", currentAtomType, &p, &n, &scatterLength) == 4)
        {
            neutronStructureFactor nsf = {p, n, scatterLength, currentAtomType};
            neutronScatteringFactors_.push_back(nsf);
        }
    }
    fclose(scatteringFile);
};

double ComputeScattering::computeS0()
{
    //significant adding errors if Debeye is not a double
    double DebeyeS0 = 0;
    for (size_t i = 0; i <= scatteringLength_.size(); ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            DebeyeS0 += scatteringLength_[i] * scatteringLength_[j];
        }
    }
    return DebeyeS0;
};

double ComputeScattering::compute_scattering(t_pbc *pbc, double q, Selection sel)
{
    //significant adding errors if Debeye is not a double
    double Debeye = 0;
    const int posCount = sel.posCount();
    for (int i = 0; i < posCount; ++i)
    {
        const SelectionPosition &pi = sel.position(i);
        auto index_i = pi.refId();
        auto scatterLength_i = getScatterLength(index_i);
        for (int j = 0; j < i; ++j)
        {
            const SelectionPosition &pj = sel.position(j);
            auto index_j = pj.refId();
            auto scatterLength_j = getScatterLength(index_j);
            double  formFactor = scatterLength_i * scatterLength_j;
            rvec  dx;
            if (pbc != nullptr)
            {
                pbc_dx(pbc, pi.x(), pj.x(), dx);
            }
            else
            {
                rvec_sub(pi.x(), pj.x(), dx);
            }
            double qDist = norm(dx) * q;
            double sinQDist = sin(qDist);
            Debeye += formFactor * (sinQDist / qDist);
        }
    }
    return Debeye;
};

real ComputeScattering::getScatterLength(int i)
{
    return scatteringLength_[atomIndex_[i]];
};

ComputeScattering makeScattering(const Selection &sel, const TopologyInformation &top)
{
    //const ArrayRef<const int> idx=sel.atomIndices();
    std::vector<int> atomIndex;
    std::vector<real> scatteringLength;
    std::vector<neutronStructureFactor> neutronScatteringFactors_;
    readScatteringFactors(neutronScatteringFactors_);
    const t_atoms &atoms = top.topology()->atoms;
    for (int i = 0; i < sel.posCount(); ++i)
    {
        const SelectionPosition &pi = sel.position(i);
        const t_atom atom = atoms.atom[i];
        for (auto scatter : neutronScatteringFactors_)
        {
            if (scatter.element == atom.elem)
            {
                scatteringLength.push_back(scatter.scatterLength);
            }
        }
        atomIndex.push_back(pi.refId());

    }
    if (atomIndex.size() != scatteringLength.size())
    {
        GMX_THROW(gmx::InvalidInputError("FAIL\n"));
    }
    return ComputeScattering(atomIndex, scatteringLength);
};

}; //namespace gmx
