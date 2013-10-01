/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifndef _waxs_debye_force_c_h
#define _waxs_debye_force_c_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "gromacs/utility/real.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/pbcutil//pbc.h"
#include "gromacs/legacyheaders/types/mdatom.h"
#include "gromacs/legacyheaders/types/forcerec.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Function that reads WAXS data from text files
 * and sets up parameters.
 *
 * \param[in]  fplog      File pointer for log file information
 * \param[in]  sfactor    Structure factors
 * \param[in]  waxs_ref   Reference S(q)
 * \param[in]  waxs_diff  Difference S(q), may be NULL
 * \param[out] waxs_out   Plot file for Scalc(q,t) (can be large)
 * \param[out] waxs_alpha Plot file for alpha(t)
 * \param[in]  cr         Communication data
 * \param[in]  ir         Parameters for algorithm
 *
 * \return the waxs_debye_force structure.
 */
waxs_debye_force_t waxs_debye_force_init(FILE             *fplog,
                                         const char       *sfactor,
                                         const char       *waxs_ref,
                                         const char       *waxs_diff,
                                         const char       *waxs_out,
                                         const char       *waxs_alpha,
                                         const t_commrec  *cr,
                                         const t_inputrec *ir);
/*! \brief
 * Clean up the contents of the structure
 */
void waxs_debye_force_done(waxs_debye_force_t wdf);

/*! \brief
 * Update the alpha parameter, either by integration or
 * using the bisection algorithm described in the manual.
 */
void waxs_debye_force_update_alpha(waxs_debye_force_t wdf);

/*! \brief
 * Evaluates scattering, compares to experimental data
 * and applies forces.
 *
 * \param[in]  fplog       File pointer for log file information
 * \param[in]  wdf         The data structure containing all the stuff
 * \param[in]  nbonds      Number of WAXS pairs
 * \param[in]  iatoms      Array of interactions
 * \param[in]  forceparams Array of force parameters
 * \param[in]  x           Coordinates
 * \param[out] f           Forces
 * \param[in]  pbc         Periodicity information (IS THIS NEEDED?)
 * \param[in]  cr          Communication data
 * \param[out] nrnb        Flop accounting
 *
 * \return The WAXS-Debye energy
 */
real waxs_debye_force_calc(FILE              *fplog,
                           waxs_debye_force_t wdf,
                           const int          nbonds,
                           const t_iatom      iatoms[],
                           const t_iparams    forceparams[],
                           rvec               x[],
                           rvec               f[],
                           const t_pbc       *pbc,
                           const t_commrec   *cr,
                           t_nrnb            *nrnb);

#ifdef __cplusplus
}
#endif

#endif /* _waxs_debye_force_c_h check */
