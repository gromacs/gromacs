/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Implements amplitude lookup for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfittingamplitudelookup.h"

#include <algorithm>
#include <iterator>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
/********************************************************************
 * DensityFittingAmplitudeLookup::Impl
 */

class DensityFittingAmplitudeLookupImpl
{
public:
    DensityFittingAmplitudeLookupImpl()                                         = default;
    DensityFittingAmplitudeLookupImpl(const DensityFittingAmplitudeLookupImpl&) = default;
    virtual ~DensityFittingAmplitudeLookupImpl()                                = default;

    virtual const std::vector<real>& operator()(ArrayRef<const real> /*chargeA*/,
                                                ArrayRef<const real> /*massT*/,
                                                ArrayRef<const int> localIndex) = 0;
    virtual std::unique_ptr<DensityFittingAmplitudeLookupImpl> clone()          = 0;
};

namespace
{
class UnitAmplitudes final : public DensityFittingAmplitudeLookupImpl
{
public:
    UnitAmplitudes()                      = default;
    UnitAmplitudes(const UnitAmplitudes&) = default;
    ~UnitAmplitudes() override            = default;
    std::unique_ptr<DensityFittingAmplitudeLookupImpl> clone() override;
    const std::vector<real>&                           operator()(ArrayRef<const real> /*chargeA*/,
                                        ArrayRef<const real> /*massT*/,
                                        ArrayRef<const int> localIndex) override;

private:
    std::vector<real> amplitude_;
};

std::unique_ptr<DensityFittingAmplitudeLookupImpl> UnitAmplitudes::clone()
{
    return std::make_unique<UnitAmplitudes>(*this);
};

const std::vector<real>& UnitAmplitudes::operator()(ArrayRef<const real> /*chargeA*/,
                                                    ArrayRef<const real> /*massT*/,
                                                    ArrayRef<const int> localIndex)
{
    if (amplitude_.size() != localIndex.size())
    {
        amplitude_.resize(localIndex.size(), 1.);
    }

    return amplitude_;
}

class ChargesAsAmplitudes final : public DensityFittingAmplitudeLookupImpl
{
public:
    ChargesAsAmplitudes()                           = default;
    ChargesAsAmplitudes(const ChargesAsAmplitudes&) = default;
    ~ChargesAsAmplitudes() override                 = default;
    std::unique_ptr<DensityFittingAmplitudeLookupImpl> clone() override;
    const std::vector<real>&                           operator()(ArrayRef<const real> chargeA,
                                        ArrayRef<const real> /*massT*/,
                                        ArrayRef<const int> localIndex) override;

private:
    std::vector<real> amplitude_;
};

std::unique_ptr<DensityFittingAmplitudeLookupImpl> ChargesAsAmplitudes::clone()
{
    return std::make_unique<ChargesAsAmplitudes>(*this);
};

const std::vector<real>& ChargesAsAmplitudes::operator()(ArrayRef<const real> chargeA,
                                                         ArrayRef<const real> /*massT*/,
                                                         ArrayRef<const int> localIndex)
{
    if (amplitude_.size() != localIndex.size())
    {
        amplitude_.resize(localIndex.size());
    }

    std::transform(std::begin(localIndex),
                   std::end(localIndex),
                   std::begin(amplitude_),
                   [&chargeA](gmx::Index index) { return chargeA[index]; });
    return amplitude_;
}

class MassesAsAmplitudes final : public DensityFittingAmplitudeLookupImpl
{
public:
    MassesAsAmplitudes()                          = default;
    MassesAsAmplitudes(const MassesAsAmplitudes&) = default;
    ~MassesAsAmplitudes() override                = default;
    std::unique_ptr<DensityFittingAmplitudeLookupImpl> clone() override;
    const std::vector<real>&                           operator()(ArrayRef<const real> /*chargeA*/,
                                        ArrayRef<const real> massT,
                                        ArrayRef<const int>  localIndex) override;

private:
    std::vector<real> amplitude_;
};

std::unique_ptr<DensityFittingAmplitudeLookupImpl> MassesAsAmplitudes::clone()
{
    return std::make_unique<MassesAsAmplitudes>(*this);
};

const std::vector<real>& MassesAsAmplitudes::operator()(ArrayRef<const real> /*chargeA*/,
                                                        ArrayRef<const real> massT,
                                                        ArrayRef<const int>  localIndex)
{
    if (amplitude_.size() != localIndex.size())
    {
        amplitude_.resize(localIndex.size());
    }

    std::transform(std::begin(localIndex),
                   std::end(localIndex),
                   std::begin(amplitude_),
                   [&massT](gmx::Index index) { return massT[index]; });
    return amplitude_;
}

} // namespace

/********************************************************************
 * DensityFittingAmplitudeLookup
 */


DensityFittingAmplitudeLookup::DensityFittingAmplitudeLookup(const DensityFittingAmplitudeMethod& method)
{
    switch (method)
    {
        case DensityFittingAmplitudeMethod::Unity:
            impl_ = std::make_unique<UnitAmplitudes>();
            break;
        case DensityFittingAmplitudeMethod::Mass:
            impl_ = std::make_unique<MassesAsAmplitudes>();
            break;
        case DensityFittingAmplitudeMethod::Charge:
            impl_ = std::make_unique<ChargesAsAmplitudes>();
            break;
        default: break;
    }
}

const std::vector<real>& DensityFittingAmplitudeLookup::operator()(ArrayRef<const real> chargeA,
                                                                   ArrayRef<const real> massT,
                                                                   ArrayRef<const int>  localIndex)
{
    return (*impl_)(chargeA, massT, localIndex);
}

DensityFittingAmplitudeLookup::~DensityFittingAmplitudeLookup() = default;

DensityFittingAmplitudeLookup::DensityFittingAmplitudeLookup(const DensityFittingAmplitudeLookup& other) :
    impl_(other.impl_->clone())
{
}

DensityFittingAmplitudeLookup& DensityFittingAmplitudeLookup::operator=(const DensityFittingAmplitudeLookup& other)
{
    impl_ = other.impl_->clone();
    return *this;
}

DensityFittingAmplitudeLookup::DensityFittingAmplitudeLookup(DensityFittingAmplitudeLookup&&) noexcept = default;

DensityFittingAmplitudeLookup&
DensityFittingAmplitudeLookup::operator=(DensityFittingAmplitudeLookup&&) noexcept = default;

} // namespace gmx
