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
 * Implementation of base class for Debye Scattering
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "scattering-debye.h"

#include <cmath>
#include <cstddef>

#include <algorithm>
#include <iterator>
#include <numeric>

#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

void ComputeDebyeScattering::initPairDistHist()
{
    // Calculate possible max Index
    maxHIndex_ = static_cast<size_t>(std::floor(maxDist_ / binWidth_) + 1);
    if (sfDepenOnQ_)
    {
        sfQDependDistValues_.resize(qValues_.size());
        for (size_t i = 0; i < qValues_.size(); ++i)
        {
            sfQDependDistValues_[i].resize(maxHIndex_);
        }
    }
    else
    {
        sfDistValues_.resize(maxHIndex_);
    }
    histRValues_.resize(maxHIndex_);
    for (size_t i = 0; i < histRValues_.size(); ++i)
    {
        histRValues_[i] = binWidth_ * i + binWidth_ * 0.5;
    }
}

void ComputeDebyeScattering::getMaxDist(matrix box)
{
    RVec boxX{ box[XX] };
    RVec boxY{ box[YY] };
    RVec boxZ{ box[ZZ] };
    RVec boxDiag = (boxX + boxY + boxZ);
    maxDist_     = boxDiag.norm();
}

double ComputeDebyeScattering::getFormFactor(int i, int j, double q)
{
    return getScatteringLength(i, q) * getScatteringLength(j, q);
}

void ComputeDebyeScattering::computeIntensity()
{
    intensity_.resize(qValues_.size());

    for (size_t qi = 0; qi < qValues_.size(); ++qi)
    {
        if (qValues_[qi] == 0.0)
        {
            intensity_[qi] = computeIntensityZeroQ();
        }
        else
        {
            for (size_t hi = 0; hi < maxHIndex_; ++hi)
            {
                double qDist = qValues_[qi] * histRValues_[hi];
                if (sfDepenOnQ_)
                {
                    intensity_[qi] += sfQDependDistValues_[qi][hi] * sin(qDist) / qDist;
                }
                else
                {
                    intensity_[qi] += sfDistValues_[hi] * sin(qDist) / qDist;
                }
            }
        }
    }
}

double ComputeDebyeScattering::computeIntensityZeroQ()
{
    if (sfDepenOnQ_)
    {
        return std::accumulate(sfQDependDistValues_[0].begin(), sfQDependDistValues_[0].end(), 0.0);
    }
    else
    {
        return std::accumulate(sfDistValues_.begin(), sfDistValues_.end(), 0.0);
    }
}

void ComputeDebyeScattering::computeDirectPairDistancesHistogram(t_pbc* pbc, Selection sel)
{
    const int posCount = sel.posCount();
    for (auto i = 0; i < posCount; ++i)
    {
        const SelectionPosition& pos_i = sel.position(i);
        for (auto j = i + 1; j < posCount; ++j)
        {
            PairDistValue            pairDist{};
            const SelectionPosition& pos_j = sel.position(j);
            pairDist.atomI                 = sel.position(i).atomIndices()[0];
            pairDist.atomJ                 = sel.position(j).atomIndices()[0];
            RVec distance_ij;
            if (pbc != nullptr)
            {
                pbc_dx(pbc, pos_i.x(), pos_j.x(), distance_ij);
            }
            else
            {
                rvec_sub(pos_i.x(), pos_j.x(), distance_ij);
            }
            pairDist.distanceIJ = distance_ij.norm();
            addPairToHist(pairDist);
        }
    }
}

void ComputeDebyeScattering::computeMonteCarloPairDistancesHistogram(t_pbc*    pbc,
                                                                     Selection sel,
                                                                     float     coverage,
                                                                     int       seed)
{
    const size_t                   posCount = sel.posCount();
    DefaultRandomEngine            rng(seed);
    UniformIntDistribution<size_t> distribution(0, posCount - 1);
    auto numPairs = static_cast<size_t>(coverage * posCount * (posCount - 1) * 0.5);
    for (size_t pair = 0; pair < numPairs; ++pair)
    {
        PairDistValue pairDist{};
        size_t        rand_i = distribution(rng);
        size_t        rand_j = distribution(rng);
        if (rand_i != rand_j)
        {
            const SelectionPosition& pos_i = sel.position(rand_i);
            const SelectionPosition& pos_j = sel.position(rand_j);
            pairDist.atomI                 = sel.position(rand_i).atomIndices()[0];
            pairDist.atomJ                 = sel.position(rand_j).atomIndices()[0];
            RVec distance_ij;
            if (pbc != nullptr)
            {
                pbc_dx(pbc, pos_i.x(), pos_j.x(), distance_ij);
            }
            else
            {
                rvec_sub(pos_i.x(), pos_j.x(), distance_ij);
            }
            pairDist.distanceIJ = distance_ij.norm();
            addPairToHist(pairDist);
        }
    }
}

void ComputeDebyeScattering::clearHist()
{
    if (sfDepenOnQ_)
    {
        for (auto QDList : sfQDependDistValues_)
        {
            std::fill(QDList.begin(), QDList.end(), 0);
        }
    }
    else
    {
        std::fill(sfDistValues_.begin(), sfDistValues_.end(), 0);
    }
}

void ComputeDebyeScattering::addPairToHist(PairDistValue pair)
{
    size_t hidx_ = std::floor(pair.distanceIJ / binWidth_);
    if (sfDepenOnQ_)
    {
        for (size_t i = 0; i != qValues_.size(); ++i)
        {
            sfQDependDistValues_[i][hidx_] += getFormFactor(pair.atomI, pair.atomJ, qValues_[i]);
        }
    }
    else
    {
        sfDistValues_[hidx_] += getFormFactor(pair.atomI, pair.atomJ, 0);
    }
}

void ComputeDebyeScattering::setBinWidth(double binWidth)
{
    binWidth_ = binWidth;
}

double ComputeDebyeScattering::getIntensity(size_t qi)
{
    return intensity_[qi];
}

void ComputeDebyeScattering::addQList(std::vector<double> qList)
{
    std::copy(qList.begin(), qList.end(), std::back_inserter(qValues_));
}

}; // namespace gmx
