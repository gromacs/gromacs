/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares gmx::EnergyFrame
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ENERGYFRAME_H
#define GMX_ENERGYANALYSIS_ENERGYFRAME_H

#include <string>
#include <vector>

#include "gromacs/legacyheaders/typedefs.h"

namespace gmx
{
/*! \brief
 * Class describing an energy frame
 */
class EnergyFrame
{
    private:
        //! The time in the simulation
        double      t_;
        //! Thw step in the simulation
        gmx_int64_t step_;
        //! The energy at this time point
        double      e_;
        //! The number of statistics points
        int         nsum_;
        //! The sum of energies over nsum_ previous steps
        double      esum_;
        //! The variance over nsum_ previous steps
        double      evar_;
    public:
        /*! \brief
         * Constructor
         * \param[in] t The time
         * \param[in] step The MD step
         * \param[in] e The instantaneous energy
         * \param[in] nsum The number of statistics points
         * \param[in] esum The sum of energies over nsum_ previous steps
         * \param[in] evar The variance over nsum_ previous steps
         */
        EnergyFrame(double t, gmx_int64_t step, double e, int nsum,
                    double esum, double evar)
        {
            t_    = t;
            step_ = step;
            e_    = e;
            nsum_ = nsum;
            esum_ = esum;
            evar_ = evar;
        }
        //! Destructor
        ~EnergyFrame() {}
        //! \return the time
        double getT() { return t_; }
        //! \return the step
        gmx_int64_t getStep() { return step_; }
        //! \return the instantaneous energy
        double getE() { return e_; }
        //! \return the number of statistics points
        int getNsum() { return nsum_; }
        //! \return the sum of energies
        double getEsum() { return esum_; }
        //! \return the variance of energies
        double getEvar() { return evar_; }
};
//! Typedef for looping over EnergyFrame
typedef std::vector<EnergyFrame>::iterator EnergyFrameIterator;

}
#endif
