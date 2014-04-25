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
/*! \file
 * \brief
 * Declares gmx::EnergyAnalysis
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ENERGYTERM_H
#define GMX_ENERGYANALYSIS_ENERGYTERM_H

#include <string>
#include <vector>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/fileio/enxio.h"
#include "energyframe.h"

namespace gmx
{
/*! \brief
 * Class describing an energy term.
 */
class EnergyTerm
{
    private:
        //! Name of the energy term
        std::string              eTerm_;
        //! Unit of this energy
        std::string              eUnit_;
        //! Number of energy terms summed so far
        gmx_int64_t              nesum_;
        //! First MD step in the analysis
        gmx_int64_t              step0_;
        //! Last MD step in the analysis
        gmx_int64_t              step1_;
        //! Index in the energy array in the energy file
        unsigned int             findex_;
        //! Best estimate of the average energy so far
        double                   energy_;
        //! Best estimae of the standard deviation so far
        double                   stddev_;
        //! Start time of the analysis
        double                   t0_;
        //! End time of the analysis
        double                   t1_;
        //! Boolean indicating whether we are storing data in the vectors below
        bool                     bStoreData_;
        //! Array of energy frames
        std::vector<EnergyFrame> ef_;
        //! Total sum of energy
        double                   esumTot_;
        //! Total variance of energy
        double                   evarTot_;
        /*! Do we have exact data from the simulations?
         * (that is, did we use nstcalcenergy == 1)?
         */
        bool                bExactData_;
        //! Is the present energy term really an energy?
        bool                bIsEner_;
        //! Has the drift been computed using linear regression?
        bool                bDrift_;
        //! Linear regression parameter a
        real                a_;
        //! Linear regression parameter b
        real                b_;
        //! Correlation coefficient from regression
        real                R_;
        //! Chi-squared from regression
        real                chi2_;

    public:
        /*! \brief
         * Constructor
         * \param[in] findex File index (in the energies stored)
         * \param[in] eTerm  String describing the energy
         * \param[in] eUnit  String describing the energy unit
         */
        EnergyTerm(unsigned int findex,
                   std::string  eTerm,
                   std::string  eUnit);

        //! Destructor
        ~EnergyTerm() {};

        //! Set all variables to zero
        void reset();

        //! Return the file index to the function type stored here
        unsigned int getIndex() { return findex_; }

        //! Return the name corresponding to the energy term
        std::string getEterm() { return eTerm_; }

        //! Return the name corresponding to the energy unit
        std::string getUnit() { return eUnit_; }

        /*! \brief
         * Tell the class to store or not to store data
         * \param[in] bStoreData Boolean
         */
        void setStoreData(bool bStoreData) { bStoreData_ = bStoreData; }

        //! Return the store data variable
        bool getStoreData() { return bStoreData_; }

        //! Is this a true energy or e.g. Temperature
        bool isEner() { return bIsEner_; }

        //! Is the exact data available?
        bool bExactData() { return bExactData_; }

        //! \return iterator to begin looping over energy frames
        EnergyFrameIterator beginEF() { return ef_.begin(); }

        //! \return iterator to end looping over energy frames
        EnergyFrameIterator endEF() { return ef_.end(); }

        /*! \brief
         * Search an energy frame corresponding to a certain step
         * \param[in] iframe The index in the array
         * \return the actual EnergyFrameIterator, or endEF() if not found
         */
        EnergyFrameIterator searchEF(gmx_int64_t iframe);

        /*! \brief
         * Add data to this Energy Term
         * \param[in] t    The time in the simulation
         * \param[in] step The simulation step
         * \param[in] nsum The number of intermediate steps for the sums
         * \param[in] esum The sum of energies over the last nsum steps
         * \param[in] evar The variance of the energies over the last nsum steps
         * \param[in] e    The energy at this point in time (trajectory)
         */
        void addData(double t, gmx_int64_t step, int nsum,
                     double esum, double evar, double e);

        //! Return the average energy
        double average();

        //! Return the standard deviation
        double standardDeviation();

        /*! \brief
         * Return an error estimate based on block averaging.
         * Requires that the energies have been stored.
         * \param[in] nb Number of blocks
         * \return the error estimate
         */
        double errorEstimate(unsigned int nb);

        /*! \brief
         * Calculate the drift - can only be done when the data is stored.
         * This is done by fitting the data to a line y = ax + b.
         * \return whether the drift was computed or not
         */
        bool calculateDrift();

        /*! \brief
         * Remove the drift - can only be done when the data is stored.
         * If necessary the drift is calculated first.
         * \return whether the drift was removed or not
         */
        bool removeDrift();

        //! Return the slope of the linear regression
        double driftA() { return a_; }

        //! Return the intercept of the linear regression
        double driftB() { return b_; }

        //! Return the correlation coefficient (R) of the linear regression
        double driftR() { return R_; }

        //! Return the chi-squared of the linear regression
        double driftChi2() { return chi2_; }

        //! Return the number of points stored
        gmx_int64_t nEnergy() { return ef_.size(); }

        //! Return the length of the data set in time
        double timeSpan() { return t1_ - t0_; }

        //! Return the begin time of the data set
        double timeBegin() { return t0_; }

        //! Return the end time of the data set
        double timeEnd() { return t1_; }

        //! Return the length of the data set in time
        gmx_int64_t nSteps() { return 1 + (step1_ - step0_); }

        //! Return the begin time of the data set
        gmx_int64_t stepBegin() { return step0_; }

        //! Return the end time of the data set
        gmx_int64_t stepEnd() { return step1_; }
};

//! Typedef for looping over EnergyTerm
typedef std::vector<EnergyTerm>::iterator EnergyTermIterator;

}
#endif
