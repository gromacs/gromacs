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
#ifndef _waxs_debye_force_h
#define _waxs_debye_force_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <math.h>
#include <boost/scoped_ptr.hpp>
#include "gromacs/utility/real.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/file.h"

namespace gmx
{
/*! \libinternal \brief
 * Implements bias force towards experimental WAXS data using the
 * Debye equation.
 *
 * \author Alexander Bjorling <alexander.bjorling@chem.gu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
class WaxsDebyeForce
{
    private:
        //! Scattering vectors, read from file
        std::vector<double>    q_;

        //! Reference data, read from file
        std::vector<double>    Sref_;

        //! Difference data, optionally read from file
        std::vector<double>    dSexp_;

        //! Computed scattering function in simulation
        std::vector<double>    Scalc_;

        //! Vector of experimental deviations, sigma(q) = Ierr(q) / I(q) + 1
        std::vector<double>    sigma_;

        //! Internal help vector
        std::vector<double>    A_;

        //! 2 D Array with scattering factors per atom type, and per q,
        //! corresponding to the q array read from the reference data.
        //! This is read from sfactor.dat (or other file if specified).
        std::vector<std::vector<double> > F_;

        //! Scattering factor at q = 0
        std::vector<double> F0_;

        //! Boundaries in real space for computing Scalc
        double                 rmin2_, rmax2_;

        //! Force constant
        double                 kWaxs_;

        //! Boundaries for the occupancy (see manual)
        double                 alpha_min_, alpha_max_;

        //! Mode of optimizing alpha (see manual)
        int                    alpha_mode_;

        //! Actual alpha, change in alpha and step size
        double                 alpha_, alpha_dot_, alpha_step_;

        //! Force on alpha and "mass" of alpha variable for alpha dynamics
        double                 Falpha_, Malpha_;

        //! Total WAXS energy as computed last time around
        double                 vtot_;

        //! Step counter and output control variable
        int                    step_, nstout_, nstcalc_;

        //! Pointers to gmx::File variables for writing output.
        //! Closed in destructor if open.
        boost::scoped_ptr<File> fp_out_, fp_alpha_;

        /*! \brief
         * Internal function for reading data files.
         *
         * \param[in]  waxs_ref  Reference S(q)
         * \param[in]  waxs_diff Difference S(q)
         */
        void readData(const char *waxs_ref,
                      const char *waxs_diff,
                      FILE       *fplog);


    public:
        /*! \brief
         * Constructor that reads WAXS data from text files
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
         */
        WaxsDebyeForce(FILE             *fplog,
                       const char       *sfactor,
                       const char       *waxs_ref,
                       const char       *waxs_diff,
                       const char       *waxs_out,
                       const char       *waxs_alpha,
                       const t_commrec  *cr,
                       const t_inputrec *ir);

        //! Destructor.
        ~WaxsDebyeForce() {}

        //! Close open files if necessary
        void finalize();

        /*! \brief
         * Evaluates scattering, compares to experimental data
         * and applies forces.
         *
         * \param[in]  fplog       File pointer for log file information
         * \param[in]  nbonds      Number of WAXS pairs
         * \param[in]  iatoms      Array of interactions
         * \param[in]  forceparams Array of force parameters
         * \param[in]  x           Coordinates
         * \param[out] f           Forces
         * \param[in]  pbc         Periodicity information (IS THIS NEEDED?)
         * \param[in]  cr          Communication data
         * \param[out] nrnb        Flop accounting (may be NULL)
         *
         * \return The WAXS-Debye energy
         */
        real calc(FILE            *fplog,
                  const int        nbonds,
                  const t_iatom    iatoms[],
                  const t_iparams  forceparams[],
                  rvec             x[],
                  rvec             f[],
                  const t_pbc     *pbc,
                  const t_commrec *cr,
                  t_nrnb          *nrnb);

        /*! \brief
         * Evaluates scattering at q = 0, useful for scaling. This
         * works only if the iatoms are complete. The result is summed
         * over processors.
         *
         * \param[in]  fplog       File pointer for debug output
         * \param[in]  nbonds      Number of WAXS pairs
         * \param[in]  iatoms      Array of interactions
         * \param[in]  forceparams Array of force parameters
         * \param[in]  cr          Communication data
         *
         * \return S(0)
         */
        real calcS0(FILE            *fplog,
                    const int        nbonds,
                    const t_iatom    iatoms[],
                    const t_iparams  forceparams[],
                    const t_commrec *cr);

        /*! \brief
         * Prints the simulated Scalc(q,t) every nstout
         * steps to an xvg file.
         * Prints alpha(t) to an xvg file
         */
        void dump();

        /*! \brief
         * Update the alpha parameter, either by integration or using the golden ratio
         * algorithm as described in the manual.
         */
        void updateAlpha();

        //! Return the current value of Alpha
        real getAlpha() { return alpha_; }

        //! Modify the value of Alpha
        void setAlpha(real alpha) { alpha_ = alpha; }

};
}

#endif /* _waxs_debye_force_h check */
