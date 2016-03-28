/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
 */

#ifndef GMX_MDTYPES_AWH_PARAMS_H
#define GMX_MDTYPES_AWH_PARAMS_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/mdtypes/md_enums.h"

//! Target distribution enum.
enum {
    eawhtargetCONSTANT, eawhtargetCUTOFF, eawhtargetBOLTZMANN, eawhtargetWEIGHTHIST, eawhtargetNR
};
//! String for target distribution.
extern const char *eawhtarget_names[eawhtargetNR+1];
//! Macro for target distribution string.
#define EAWHTARGET(e)  enum_name(e, eawhtargetNR, eawhtarget_names)

//! Weight histogram growth enum.
enum {
    eawhgrowthEXP_LINEAR, eawhgrowthNONE, eawhgrowthLINEAR, eawhgrowthNR
};
//! String for weight histogram growth
extern const char *eawhgrowth_names[eawhgrowthNR+1];
//! Macro for weight histogram growth string.
#define EAWHGROWTH(e)    enum_name(e, eawhgrowthNR, eawhgrowth_names)

//! Parameters for an AWH coordinate dimension.
typedef struct awhdim_params_t {
    int    pull_coord_index;          /**< Index of the pull coordinate to bias. */
    double period;                    /**< The period of this dimension (= 0 if not periodic). */
    double diffusion;                 /**< Estimated diffusion constant in units of nm^2/ps or rad^2/ps. */
    double origin;                    /**< Start value of the interval. */
    double end;                       /**< End value of the interval. */
    int    ninterval;                 /**< The number of subintervals to split the interval into along this dimension. */
    double interval_overlap;          /**< Subinterval overlap (as a fractional value in [0, 1]). */
    double coord_value_init;          /**< The initial coordinate value. */
} awhdim_params_t;

//! Parameters for an AWH bias.
typedef struct awh_params_t {
    int                ndim;               /**< Dimension of the coordinate space. */
    awhdim_params_t   *dim_params;         /**< AWH parameters per dimension. */
    int                eTarget;            /**< Type of target distribution. */
    double             target_param;       /**< Target distribution parameter (meaning depends on eTarget). */
    int                eGrowth;            /**< How the biasing histogram grows. */
    int                bUser_data;         /**< Is there a user-defined initial PMF estimate and target estimate? */
    double             error_initial;      /**< Estimated initial free energy error. */
    gmx_bool           bShare;             /**< Share the bias across multiple simulations? */
} awh_params_t;

//! Parameters for a collection of AWH biases.
typedef struct awhbias_params_t {
    int                 nawh;                        /**< The number of AWH biases.*/
    awh_params_t       *awh_params;                  /**< AWH bias parameters.*/
    gmx_int64_t         seed;                        /**< Random seed.*/
    int                 nstsample_coord;             /**< Number of samples per coordinate sample (also used for PMF) */
    int                 nsamples_move_refvalue;      /**< Number of samples per moving the coordinate reference value. */
    int                 nsamples_update_free_energy; /**< Number of samples per free energy update. */
    gmx_bool            bConvolve_force;             /**< If true, convolved the bias force. Otherwise, apply harmonic potential force. */
} awhbias_params_t;

#endif /* GMX_MDTYPES_AWH_PARAMS_H */
