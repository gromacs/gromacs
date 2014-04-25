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
#ifndef GMX_FEE_H
#define GMX_FEE_H

#include <stdio.h>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/oenv.h"
#include "enerdata.h"
#include "energyanalysis.h"
#include "simple.h"

namespace gmx
{
class FreeEnergyEstimate : public SimpleEnergy
{
    private:
        //! Second energy file
        std::string eneFile2_;

        //! Running average for anlysis
        std::string runAverFile_;

        //! Reference temperature
        double      reftemp_;

        int         nset, *set;
        char      **leg;
        enerdata_t *edat;
        double     *time;

    public:
        //! Constructor
        FreeEnergyEstimate() { reftemp_ = 0; };

        //! Destructor
        ~FreeEnergyEstimate() {};

        /*! \brief
         * Some option
         * \param[in] reftemp Reference Temperature
         */
        void setParameters(double reftemp) { reftemp_ = reftemp; }

        //! Set the fluctuation convergence file name
        void setEneFile2(std::string eneFile2)
        {
            eneFile2_ = eneFile2;
        }

        //! Get the fluctuation convergence file name
        std::string getEneFile2() { return eneFile2_; }

        //! Set the fluctuation convergence file name
        void setRunAverFile(std::string runAverFile)
        {
            runAverFile_ = runAverFile;
        }

        //! Get the fluctuation convergence file name
        std::string getRunAverFile() { return runAverFile_; }

        //! Does the initiation of the reading of the file
        virtual bool readInit(int nre, gmx_enxnm_t enm[]);

        //! Analyse one frame and stores the results in memory
        virtual bool readFrame(t_enxframe *fr);

        //! Finalize reading
        virtual bool readFinalize();
};

};

#endif
