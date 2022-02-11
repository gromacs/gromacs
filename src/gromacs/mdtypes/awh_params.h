/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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

/*! \libinternal \file
 *
 * \brief
 * Declares AWH parameter data types.
 *
 * Besides internal use by the AWH module, the AWH parameters are needed
 * for reading the user input (mdp) file and for reading and writing the
 * parameters to the mdrun input (tpr) file.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_AWH_PARAMS_H
#define GMX_MDTYPES_AWH_PARAMS_H

#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"

struct t_inpfile;
struct t_inputrec;
struct pull_params_t;
using warninp_t = struct warninp*;

namespace gmx
{

class ISerializer;
//! Target distribution enum.
enum class AwhTargetType : int
{
    Constant,
    Cutoff,
    Boltzmann,
    LocalBoltzmann,
    Count,
    Default = Constant
};
//! String for target distribution.
const char* enumValueToString(AwhTargetType enumValue);

//! Weight histogram growth enum.
enum class AwhHistogramGrowthType : int
{
    ExponentialLinear,
    Linear,
    Count,
    Default = ExponentialLinear
};
//! String for weight histogram growth
const char* enumValueToString(AwhHistogramGrowthType enumValue);

//! AWH potential type enum.
enum class AwhPotentialType : int
{
    Convolved,
    Umbrella,
    Count,
    Default = Convolved
};
//! String for AWH potential type
const char* enumValueToString(AwhPotentialType enumValue);

//! AWH bias reaction coordinate provider
enum class AwhCoordinateProviderType : int
{
    Pull,
    FreeEnergyLambda,
    Count,
    Default = Pull
};
//! String for AWH bias reaction coordinate provider.
const char* enumValueToString(AwhCoordinateProviderType enumValue);

class AwhDimParams
{
public:
    //! Constructor from input file.
    AwhDimParams(std::vector<t_inpfile>* inp, const std::string& prefix, warninp_t wi, bool bComment);
    //! Constructor to generate from file reading.
    explicit AwhDimParams(ISerializer* serializer);

    //! Move constructor.
    AwhDimParams(AwhDimParams&&) = default;
    //! Move assignment operator.
    AwhDimParams& operator=(AwhDimParams&&) = default;
    //! Delete copy constructor.
    AwhDimParams(const AwhDimParams&) = delete;
    //! Delete copy assignment.
    AwhDimParams& operator=(const AwhDimParams&) = delete;

    //! Which module is providing the reaction coordinate.
    AwhCoordinateProviderType coordinateProvider() const { return eCoordProvider_; }
    //! Index for reaction coordinate in provider.
    int coordinateIndex() const { return coordIndex_; }
    //! Start value for interval.
    double origin() const { return origin_; }
    //! End value for interval.
    double end() const { return end_; }
    //! Period for the dimension.
    double period() const { return period_; }
    //! Set period value dependent on state.
    void setPeriod(double period) { period_ = period; }
    //! Force constant for this dimension.
    double forceConstant() const { return forceConstant_; }
    //! Estimated diffusion constant.
    double diffusion() const { return diffusion_; }
    //! Initial value for coordinate.
    double initialCoordinate() const { return coordValueInit_; }
    //! Set initial coordinate value dependent on state.
    void setInitialCoordinate(double initialCoordinate) { coordValueInit_ = initialCoordinate; }
    //! Diameter needed to be sampled.
    double coverDiameter() const { return coverDiameter_; }
    //! Write datastructure.
    void serialize(ISerializer* serializer);

private:
    //! The module providing the reaction coordinate.
    AwhCoordinateProviderType eCoordProvider_;
    //! Index of reaction coordinate in the provider.
    int coordIndex_ = 0;
    //! Start value of the interval.
    double origin_ = 0.0;
    //! End value of the interval.
    double end_ = 0.0;
    //! The period of this dimension (= 0 if not periodic).
    double period_ = 0.0;
    //! The force constant in kJ/mol/nm^2, kJ/mol/rad^2
    double forceConstant_ = 0.0;
    //! Estimated diffusion constant in units of nm^2/ps or rad^2/ps or ps^-1.
    double diffusion_ = 0.0;
    //! The initial coordinate value.
    double coordValueInit_ = 0.0;
    //! The diameter that needs to be sampled around a point before it is considered covered.
    double coverDiameter_ = 0.0;
};

class AwhBiasParams
{
public:
    //! Constructor from input file.
    AwhBiasParams(std::vector<t_inpfile>* inp, const std::string& prefix, warninp_t wi, bool bComment);
    //! Constructor to generate from file reading.
    explicit AwhBiasParams(ISerializer* serializer);

    //! Move constructor.
    AwhBiasParams(AwhBiasParams&&) = default;
    //! Move assignment operator.
    AwhBiasParams& operator=(AwhBiasParams&&) = default;
    //! Delete copy constructor.
    AwhBiasParams(const AwhBiasParams&) = delete;
    //! Delete copy assignment.
    AwhBiasParams& operator=(const AwhBiasParams&) = delete;

    //! Which target distribution is searched.
    AwhTargetType targetDistribution() const { return eTarget_; }
    //! Beta scaling to reach target distribution.
    double targetBetaScaling() const { return targetBetaScaling_; }
    //! Cutoff for target.
    double targetCutoff() const { return targetCutoff_; }
    //! Which kind of growth to use.
    AwhHistogramGrowthType growthType() const { return eGrowth_; }
    //! User provided PMF estimate.
    bool userPMFEstimate() const { return bUserData_; }
    //! Estimated initial free energy error in kJ/mol.
    double initialErrorEstimate() const { return errorInitial_; }
    //! Dimensions of coordinate space.
    int ndim() const { return dimParams_.size(); }
    //! Number of groups to share this bias with.
    int shareGroup() const { return shareGroup_; }
    //! If the simulation starts with equilibrating histogram.
    bool equilibrateHistogram() const { return equilibrateHistogram_; }
    //! Access to dimension parameters.
    ArrayRef<AwhDimParams> dimParams() { return dimParams_; }
    //! Const access to dimension parameters.
    ArrayRef<const AwhDimParams> dimParams() const { return dimParams_; }
    //! Write datastructure.
    void serialize(ISerializer* serializer);

private:
    //! AWH parameters per dimension.
    std::vector<AwhDimParams> dimParams_;
    //! Type of target distribution.
    AwhTargetType eTarget_;
    //! Beta scaling value for Boltzmann type target distributions.
    double targetBetaScaling_;
    //! Free energy cutoff value for cutoff type target distribution in kJ/mol.
    double targetCutoff_;
    //! How the biasing histogram grows.
    AwhHistogramGrowthType eGrowth_;
    //! Is there a user-defined initial PMF estimate and target estimate?
    bool bUserData_;
    //! Estimated initial free energy error in kJ/mol.
    double errorInitial_;
    //! When >0, the bias is shared with biases of the same group and across multiple simulations when shareBiasMultisim=true
    int shareGroup_;
    //! True if the simulation starts out by equilibrating the histogram.
    bool equilibrateHistogram_;
};
/*! \internal
 * \brief Structure holding parameter information for AWH.
 */
class AwhParams
{
public:
    //! Constructor from input file.
    AwhParams(std::vector<t_inpfile>* inp, warninp_t wi);
    //! Constructor used to generate awh parameter from file reading.
    explicit AwhParams(ISerializer* serializer);

    //! Move constructor.
    AwhParams(AwhParams&&) = default;
    //! Move assignment operator.
    AwhParams& operator=(AwhParams&&) = default;
    //! Delete copy constructor.
    AwhParams(const AwhParams&) = delete;
    //! Delete copy assignment.
    AwhParams& operator=(const AwhParams&) = delete;

    //! Get number of biases.
    int numBias() const { return awhBiasParams_.size(); }
    //! Get access to bias parameters.
    ArrayRef<AwhBiasParams> awhBiasParams() { return awhBiasParams_; }
    //! Const access to bias parameters.
    ArrayRef<const AwhBiasParams> awhBiasParams() const { return awhBiasParams_; }
    //! What king of potential is being used. \todo should use actual enum class.
    AwhPotentialType potential() const { return potentialEnum_; }
    //! Seed used for starting AWH.
    int64_t seed() const { return seed_; }
    //! Output step interval.
    int nstout() const { return nstOut_; }
    //! Number of samples per coordinate sample.
    int nstSampleCoord() const { return nstSampleCoord_; }
    //! Number of samples per free energy update.
    int numSamplesUpdateFreeEnergy() const { return numSamplesUpdateFreeEnergy_; }
    //! If biases are shared in multisim.
    bool shareBiasMultisim() const { return shareBiasMultisim_; }
    //! Serialize awh parameters.
    void serialize(ISerializer* serializer);

private:
    //! AWH bias parameters.
    std::vector<AwhBiasParams> awhBiasParams_;
    //! Random seed.
    int64_t seed_;
    //! Output step interval.
    int nstOut_;
    //! Number of samples per coordinate sample (also used for PMF)
    int nstSampleCoord_;
    //! Number of samples per free energy update.
    int numSamplesUpdateFreeEnergy_;
    //! Type of potential.
    AwhPotentialType potentialEnum_;
    //! Whether to share biases with shareGroup>0 between multi-simulations.
    bool shareBiasMultisim_;
};

} // namespace gmx

#endif /* GMX_MDTYPES_AWH_PARAMS_H */
