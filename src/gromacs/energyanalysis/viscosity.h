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
#ifndef GMX_VISCOSITY_H
#define GMX_VISCOSITY_H

#include "gromacs/fileio/enxio.h"
#include "simple.h"

namespace gmx
{
class Viscosity : public SimpleEnergy
{
    private:
        //! Filename for Einstein viscosity
        std::string fnEinstein_;
        //! Filename for Einstein viscosity integral
        std::string fnEinsteinIntegral_;
        /*! \brief
         * Compute viscosity using Einstein relation, see GROMACS manual.
         * The output is written to two files.
         * \param[in] nsets Number of sum data sets
         * \param[in] sum The sum data sets
         * \param[in] V The average volume
         * \param[in] T The average temperature
         */
        void doEinstein(int nsets, real **sum,
                        real V, real T);
    public:
        //! Constructor
        Viscosity() {};

        //! Destructor must be virtual since there are other virtual functions
        virtual ~Viscosity() {};

        //! Set output file name for Einstein
        void setFnEinstein(std::string fn) { fnEinstein_ = fn; }

        std::string getFnEinstein() { return fnEinstein_; }

        //! Set output file name for Einstein integral
        void setFnEinsteinIntegral(std::string fn) { fnEinsteinIntegral_ = fn; }

        std::string getFnEinsteinIntegral() { return fnEinsteinIntegral_; }

        //! Does the initiation of the reading of the file
        virtual bool initAnalysis(int nre, gmx_enxnm_t enm[]);

        //! Analyse one frame and stores the results in memory
        virtual bool addAnalysisFrame(t_enxframe *fr);

        //! Finalize reading
        virtual bool finalizeAnalysis();
};

};

#endif
