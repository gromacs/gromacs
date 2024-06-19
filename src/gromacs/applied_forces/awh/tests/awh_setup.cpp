/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
#include "gmxpre.h"

#include "awh_setup.h"

#include <cmath>
#include <cstdint>

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/bias.h"
#include "gromacs/applied_forces/awh/correlationgrid.h"
#include "gromacs/applied_forces/awh/pointstate.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
class ISerializer;

namespace test
{

using ::testing::Eq;
using ::testing::Pointwise;

std::vector<char> awhDimParamSerialized(AwhCoordinateProviderType inputCoordinateProvider,
                                        int                       inputCoordIndex,
                                        double                    inputOrigin,
                                        double                    inputEnd,
                                        double                    inputPeriod,
                                        double                    inputDiffusion)
{
    AwhCoordinateProviderType eCoordProvider = inputCoordinateProvider;
    int                       coordIndex     = inputCoordIndex;
    double                    forceConstant  = 10;
    double                    period         = inputPeriod;
    double                    diffusion      = inputDiffusion;
    double                    origin         = inputOrigin;
    double                    end            = inputEnd;
    double                    coordValueInit = inputOrigin;
    double                    coverDiameter  = 0;

    gmx::InMemorySerializer serializer;
    serializer.doEnumAsInt(&eCoordProvider);
    serializer.doInt(&coordIndex);
    serializer.doDouble(&origin);
    serializer.doDouble(&end);
    serializer.doDouble(&period);
    serializer.doDouble(&forceConstant);
    serializer.doDouble(&diffusion);
    serializer.doDouble(&coordValueInit);
    serializer.doDouble(&coverDiameter);
    return serializer.finishAndGetBuffer();
}

/*! \internal \brief
 * Prepare a memory buffer with serialized AwhBiasParams.
 *
 * \param[in] eawhgrowth Way to grow potential.
 * \param[in] beta Value for 1/(kB*T).
 * \param[in] inputErrorScaling Factor for initial error scaling.
 * \param[in] dimensionParameterBuffers Buffers containing the dimension parameters.
 * \param[in] shareGroup share group for, potentially, sharing the bias between simulations
 * \param[in] inputUserData If there is a user provided PMF estimate.
 * \param[in] eTargetType Target distribution type.
 * \param[in] scaleTargetByMetric Whether to scale the target distribution based on the friction metric.
 */
static std::vector<char> awhBiasParamSerialized(AwhHistogramGrowthType            eawhgrowth,
                                                double                            beta,
                                                double                            inputErrorScaling,
                                                ArrayRef<const std::vector<char>> dimensionParameterBuffers,
                                                int                               shareGroup,
                                                bool                              inputUserData,
                                                AwhTargetType                     eTargetType,
                                                bool scaleTargetByMetric)
{
    int                    ndim                     = dimensionParameterBuffers.size();
    double                 targetBetaScaling        = 0;
    double                 targetCutoff             = 0;
    AwhHistogramGrowthType eGrowth                  = eawhgrowth;
    double                 growthFactor             = 3.0;
    bool                   bUserData                = inputUserData;
    double                 errorInitial             = inputErrorScaling / beta;
    bool                   equilibrateHistogram     = false;
    double                 targetMetricScalingLimit = 10;

    gmx::InMemorySerializer serializer;
    serializer.doEnumAsInt(&eTargetType);
    serializer.doDouble(&targetBetaScaling);
    serializer.doDouble(&targetCutoff);
    serializer.doEnumAsInt(&eGrowth);
    serializer.doDouble(&growthFactor);
    int temp = static_cast<int>(bUserData);
    serializer.doInt(&temp);
    serializer.doBool(&scaleTargetByMetric);
    serializer.doDouble(&targetMetricScalingLimit);
    serializer.doDouble(&errorInitial);
    serializer.doInt(&ndim);
    serializer.doInt(&shareGroup);
    serializer.doBool(&equilibrateHistogram);

    auto awhDimBuffer  = awhDimParamSerialized();
    auto awhBiasBuffer = serializer.finishAndGetBuffer();
    for (const auto& dimParamBuffer : dimensionParameterBuffers)
    {
        awhBiasBuffer.insert(awhBiasBuffer.end(), dimParamBuffer.begin(), dimParamBuffer.end());
    }
    return awhBiasBuffer;
}

/*! \internal \brief
 * Prepare a memory buffer with serialized AwhParams.
 *
 * \param[in] eawhgrowth Way to grow potential.
 * \param[in] eawhpotential Which potential to use.
 * \param[in] beta Value for 1/(kB*T).
 * \param[in] inputErrorScaling Factor for initial error scaling.
 * \param[in] inputSeed Seed value to use.
 * \param[in] dimensionParameterBuffers Buffers containing the dimension parameters.
 * \param[in] biasShareGroup share group for, potentially, sharing the bias over simulations
 * \param[in] inputUserData If there is a user provided PMF estimate.
 * \param[in] eTargetType Target distribution type.
 * \param[in] scaleTargetByMetric Whether to scale the target distribution based on the friction metric.
 */
static std::vector<char> awhParamSerialized(AwhHistogramGrowthType            eawhgrowth,
                                            AwhPotentialType                  eawhpotential,
                                            double                            beta,
                                            double                            inputErrorScaling,
                                            int64_t                           inputSeed,
                                            ArrayRef<const std::vector<char>> dimensionParameterBuffers,
                                            int                               biasShareGroup,
                                            bool                              inputUserData,
                                            AwhTargetType                     eTargetType,
                                            bool                              scaleTargetByMetric)
{
    int              numBias                    = 1;
    int64_t          seed                       = inputSeed;
    int              nstOut                     = 0;
    int              nstSampleCoord             = 1;
    int              numSamplesUpdateFreeEnergy = 10;
    AwhPotentialType ePotential                 = eawhpotential;
    bool             shareBiasMultisim          = false;

    gmx::InMemorySerializer serializer;
    serializer.doInt(&numBias);
    serializer.doInt(&nstOut);
    serializer.doInt64(&seed);
    serializer.doInt(&nstSampleCoord);
    serializer.doInt(&numSamplesUpdateFreeEnergy);
    serializer.doEnumAsInt(&ePotential);
    serializer.doBool(&shareBiasMultisim);

    auto awhParamBuffer = serializer.finishAndGetBuffer();
    auto awhBiasBuffer  = awhBiasParamSerialized(eawhgrowth,
                                                beta,
                                                inputErrorScaling,
                                                dimensionParameterBuffers,
                                                biasShareGroup,
                                                inputUserData,
                                                eTargetType,
                                                scaleTargetByMetric);

    awhParamBuffer.insert(awhParamBuffer.end(), awhBiasBuffer.begin(), awhBiasBuffer.end());

    return awhParamBuffer;
}

AwhTestParameters::AwhTestParameters(ISerializer* serializer) : awhParams(serializer, false, false)
{
}
/*! \brief
 * Helper function to set up the C-style AWH parameters for the test.
 *
 * Builds the test input data from serialized data.
 */
AwhTestParameters getAwhTestParameters(AwhHistogramGrowthType            eawhgrowth,
                                       AwhPotentialType                  eawhpotential,
                                       ArrayRef<const std::vector<char>> dimensionParameterBuffers,
                                       bool                              inputUserData,
                                       double                            beta,
                                       bool                              useAwhFep,
                                       double                            inputErrorScaling,
                                       int                               numFepLambdaStates,
                                       int                               biasShareGroup,
                                       AwhTargetType                     eTargetType,
                                       bool                              scaleTargetByMetric)
{
    double  convFactor = 1;
    double  k          = 1000;
    int64_t seed       = 93471803;

    auto                      awhParamBuffer = awhParamSerialized(eawhgrowth,
                                             eawhpotential,
                                             beta,
                                             inputErrorScaling,
                                             seed,
                                             dimensionParameterBuffers,
                                             biasShareGroup,
                                             inputUserData,
                                             eTargetType,
                                             scaleTargetByMetric);
    gmx::InMemoryDeserializer deserializer(awhParamBuffer, false);
    AwhTestParameters         params(&deserializer);

    params.beta = beta;

    if (useAwhFep)
    {
        params.dimParams.emplace_back(DimParams::fepLambdaDimParams(numFepLambdaStates, params.beta));
    }
    else
    {
        params.dimParams.emplace_back(DimParams::pullDimParams(convFactor, k, params.beta));
    }
    return params;
}

TEST(SerializationTest, CanSerializeDimParams)
{
    auto                      awhDimBuffer = awhDimParamSerialized();
    gmx::InMemoryDeserializer deserializer(awhDimBuffer, false);
    AwhDimParams              awhDimParams(&deserializer);
    EXPECT_EQ(awhDimParams.coordinateProvider(), AwhCoordinateProviderType::Pull);
    EXPECT_EQ(awhDimParams.coordinateIndex(), 0);
    EXPECT_FLOAT_EQ(awhDimParams.forceConstant(), 10);
    EXPECT_FLOAT_EQ(awhDimParams.period(), 0);
    EXPECT_FLOAT_EQ(awhDimParams.diffusion(), 0.34690997);
    EXPECT_FLOAT_EQ(awhDimParams.origin(), 0.5);
    EXPECT_FLOAT_EQ(awhDimParams.end(), 1.5);
    EXPECT_FLOAT_EQ(awhDimParams.initialCoordinate(), awhDimParams.origin());
    EXPECT_FLOAT_EQ(awhDimParams.coverDiameter(), 0);

    gmx::InMemorySerializer serializer;
    awhDimParams.serialize(&serializer);
    EXPECT_THAT(awhDimBuffer, Pointwise(Eq(), serializer.finishAndGetBuffer()));
}

TEST(SerializationTest, CanSerializeBiasParams)
{
    auto awhDimBuffer   = awhDimParamSerialized();
    auto awhDimArrayRef = gmx::arrayRefFromArray(&awhDimBuffer, 1);
    auto awhBiasBuffer  = awhBiasParamSerialized(
            AwhHistogramGrowthType::ExponentialLinear, 0.4, 0.5, awhDimArrayRef, 0, false, AwhTargetType::Constant, false);
    gmx::InMemoryDeserializer deserializer(awhBiasBuffer, false);
    AwhBiasParams             awhBiasParams(&deserializer, false, false);
    EXPECT_EQ(awhBiasParams.ndim(), 1);
    EXPECT_EQ(awhBiasParams.targetDistribution(), AwhTargetType::Constant);
    EXPECT_FLOAT_EQ(awhBiasParams.targetBetaScaling(), 0);
    EXPECT_FLOAT_EQ(awhBiasParams.targetCutoff(), 0);
    EXPECT_EQ(awhBiasParams.growthType(), AwhHistogramGrowthType::ExponentialLinear);
    EXPECT_EQ(awhBiasParams.scaleTargetByMetric(), false);
    EXPECT_FLOAT_EQ(awhBiasParams.targetMetricScalingLimit(), 10);
    EXPECT_EQ(awhBiasParams.userPMFEstimate(), 0);
    EXPECT_FLOAT_EQ(awhBiasParams.initialErrorEstimate(), 0.5 / 0.4);
    EXPECT_EQ(awhBiasParams.shareGroup(), 0);
    EXPECT_EQ(awhBiasParams.equilibrateHistogram(), false);
    const auto& awhDimParams = awhBiasParams.dimParams(0);
    EXPECT_EQ(awhDimParams.coordinateProvider(), AwhCoordinateProviderType::Pull);
    EXPECT_EQ(awhDimParams.coordinateIndex(), 0);
    EXPECT_FLOAT_EQ(awhDimParams.forceConstant(), 10);
    EXPECT_FLOAT_EQ(awhDimParams.period(), 0);
    EXPECT_FLOAT_EQ(awhDimParams.diffusion(), 0.34690997);
    EXPECT_FLOAT_EQ(awhDimParams.origin(), 0.5);
    EXPECT_FLOAT_EQ(awhDimParams.end(), 1.5);
    EXPECT_FLOAT_EQ(awhDimParams.initialCoordinate(), awhDimParams.origin());
    EXPECT_FLOAT_EQ(awhDimParams.coverDiameter(), 0);

    gmx::InMemorySerializer serializer;
    awhBiasParams.serialize(&serializer);
    EXPECT_THAT(awhBiasBuffer, Pointwise(Eq(), serializer.finishAndGetBuffer()));
}

TEST(SerializationTest, CanSerializeAwhParams)
{
    auto awhDimBuffer   = awhDimParamSerialized();
    auto awhDimArrayRef = gmx::arrayRefFromArray(&awhDimBuffer, 1);
    auto awhParamBuffer = awhParamSerialized(AwhHistogramGrowthType::ExponentialLinear,
                                             AwhPotentialType::Convolved,
                                             0.4,
                                             0.5,
                                             1337,
                                             awhDimArrayRef,
                                             0,
                                             false,
                                             AwhTargetType::Constant,
                                             false);
    gmx::InMemoryDeserializer deserializer(awhParamBuffer, false);
    AwhParams                 awhParams(&deserializer, false, false);
    EXPECT_EQ(awhParams.numBias(), 1);
    EXPECT_EQ(awhParams.seed(), 1337);
    EXPECT_EQ(awhParams.nstout(), 0);
    EXPECT_EQ(awhParams.nstSampleCoord(), 1);
    EXPECT_EQ(awhParams.numSamplesUpdateFreeEnergy(), 10);
    EXPECT_EQ(awhParams.potential(), AwhPotentialType::Convolved);
    EXPECT_EQ(awhParams.shareBiasMultisim(), false);
    const auto& awhBiasParams = awhParams.awhBiasParams(0);
    EXPECT_EQ(awhBiasParams.ndim(), 1);
    EXPECT_EQ(awhBiasParams.targetDistribution(), AwhTargetType::Constant);
    EXPECT_FLOAT_EQ(awhBiasParams.targetBetaScaling(), 0);
    EXPECT_FLOAT_EQ(awhBiasParams.targetCutoff(), 0);
    EXPECT_EQ(awhBiasParams.growthType(), AwhHistogramGrowthType::ExponentialLinear);
    EXPECT_EQ(awhBiasParams.userPMFEstimate(), 0);
    EXPECT_FLOAT_EQ(awhBiasParams.initialErrorEstimate(), 0.5 / 0.4);
    EXPECT_EQ(awhBiasParams.shareGroup(), 0);
    EXPECT_EQ(awhBiasParams.equilibrateHistogram(), false);
    const auto& awhDimParams = awhBiasParams.dimParams(0);
    EXPECT_EQ(awhDimParams.coordinateProvider(), AwhCoordinateProviderType::Pull);
    EXPECT_EQ(awhDimParams.coordinateIndex(), 0);
    EXPECT_FLOAT_EQ(awhDimParams.forceConstant(), 10);
    EXPECT_FLOAT_EQ(awhDimParams.period(), 0);
    EXPECT_FLOAT_EQ(awhDimParams.diffusion(), 0.34690997);
    EXPECT_FLOAT_EQ(awhDimParams.origin(), 0.5);
    EXPECT_FLOAT_EQ(awhDimParams.end(), 1.5);
    EXPECT_FLOAT_EQ(awhDimParams.initialCoordinate(), awhDimParams.origin());
    EXPECT_FLOAT_EQ(awhDimParams.coverDiameter(), 0);

    gmx::InMemorySerializer serializer;
    awhParams.serialize(&serializer);
    EXPECT_THAT(awhParamBuffer, Pointwise(Eq(), serializer.finishAndGetBuffer()));
}

} // namespace test
} // namespace gmx
