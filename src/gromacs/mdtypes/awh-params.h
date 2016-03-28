/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

//! Target distribution enum.
enum {
    eawhtargetCONSTANT, eawhtargetCUTOFF, eawhtargetBOLTZMANN, eawhtargetLOCALBOLTZMANN, eawhtargetNR
};
//! String for target distribution.
extern const char *eawhtarget_names[eawhtargetNR+1];
//! Macro for target distribution string.
#define EAWHTARGET(e)  enum_name(e, gmx::eawhtargetNR, gmx::eawhtarget_names)

//! Weight histogram growth enum.
enum {
    eawhgrowthEXP_LINEAR, eawhgrowthLINEAR, eawhgrowthNR
};
//! String for weight histogram growth
extern const char *eawhgrowth_names[eawhgrowthNR+1];
//! Macro for weight histogram growth string.
#define EAWHGROWTH(e)    enum_name(e, gmx::eawhgrowthNR, gmx::eawhgrowth_names)

//! AWH potential type enum.
enum {
    eawhpotentialCONVOLVED, eawhpotentialUMBRELLA, eawhpotentialNR
};
//! String for AWH potential type
extern const char *eawhpotential_names[eawhpotentialNR+1];
//! Macro for AWH potential type string.
#define EAWHPOTENTIAL(e)    enum_name(e, gmx::eawhpotentialNR, gmx::eawhpotential_names)

//! AWH bias reaction coordinate provider
enum {
    eawhcoordproviderPULL, eawhcoordproviderNR
};
//! String for AWH bias reaction coordinate provider.
extern const char *eawhcoordprovider_names[eawhcoordproviderNR+1];
//! Macro for AWH bias reaction coordinate provider.
#define EAWHCOORDPROVIDER(e)    enum_name(e, gmx::eawhcoordproviderNR, gmx::eawhcoordprovider_names)

/*! \cond INTERNAL */

//! Parameters for an AWH coordinate dimension.
struct AwhDimParams
{
    int    eCoordProvider; /**< The module providing the reaction coordinate. */
    int    coordIndex;     /**< Index of reaction coordinate in the provider. */
    double origin;         /**< Start value of the interval. */
    double end;            /**< End value of the interval. */
    double period;         /**< The period of this dimension (= 0 if not periodic). */
    double forceConstant;  /**< The force constant in kJ/mol/nm^2, kJ/mol/rad^2 */
    double diffusion;      /**< Estimated diffusion constant in units of nm^2/ps or rad^2/ps. */
    double coordValueInit; /**< The initial coordinate value. */
    double coverDiameter;  /**< The diameter that needs to be sampled around a point before it is considered covered. */
};

//! Parameters for an AWH bias.
struct AwhBiasParams
{
    // TODO: Turn dimParams into a std::vector when moved into AWH module
    int           ndim;                  /**< Dimension of the coordinate space. */
    AwhDimParams *dimParams;             /**< AWH parameters per dimension. */
    int           eTarget;               /**< Type of target distribution. */
    double        targetBetaScaling;     /**< Beta scaling value for Boltzmann type target distributions. */
    double        targetCutoff;          /**< Free energy cutoff value for cutoff type target distribution in kJ/mol.*/
    int           eGrowth;               /**< How the biasing histogram grows. */
    int           bUserData;             /**< Is there a user-defined initial PMF estimate and target estimate? */
    double        errorInitial;          /**< Estimated initial free energy error in kJ/mol. */
    int           shareGroup;            /**< When >0, the bias is shared with biases of the same group and across multiple simulations when shareBiasMultisim=true */
    gmx_bool      equilibrateHistogram;  /**< True if the simulation starts out by equilibrating the histogram.  */
};

//! Parameters for AWH.
struct AwhParams
{
    // TODO: Turn awhBiasParams into a std::vector when moved into AWH module
    int            numBias;                    /**< The number of AWH biases.*/
    AwhBiasParams *awhBiasParams;              /**< AWH bias parameters.*/
    gmx_int64_t    seed;                       /**< Random seed.*/
    int            nstOut;                     /**< Output step interval.*/
    int            nstSampleCoord;             /**< Number of samples per coordinate sample (also used for PMF) */
    int            numSamplesUpdateFreeEnergy; /**< Number of samples per free energy update. */
    int            ePotential;                 /**< Type of potential. */
    gmx_bool       shareBiasMultisim;          /**< When true, share biases with shareGroup>0 between multi-simulations */
};

/*! \endcond */

}      // namespace gmx

#endif /* GMX_MDTYPES_AWH_PARAMS_H */
