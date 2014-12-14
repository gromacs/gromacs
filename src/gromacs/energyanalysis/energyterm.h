/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
 * \brief
 * Declares gmx::energyanalysis::EnergyTerm
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ENERGYTERM_H
#define GMX_ENERGYANALYSIS_ENERGYTERM_H

#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "energyframe.h"

namespace gmx
{

namespace energyanalysis
{

//! Typedef for looping over EnergyFrame
typedef std::vector<EnergyFrame>::const_iterator EnergyFrameIterator;

/*! \libinternal
 * \brief
 * Class describing the whole time series of an energy term.
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
        //! Best estimate of the standard deviation so far
        double                   stddev_;
        //! Start time of the analysis
        double                   t0_;
        //! End time of the analysis
        double                   t1_;
        //! Boolean whether firstframe has been read
        bool                     firstFrameRead_;
        //! Boolean indicating whether we are storing data in the vectors below
        bool                     storeData_;
        //! Array of energy frames
        std::vector<EnergyFrame> ef_;
        //! Total sum of energy
        double                   esumTot_;
        //! Total variance of energy
        double                   evarTot_;
        //! Is the present energy term really an energy?
        bool                     isEner_;

    public:
        /*! \brief
         * Constructor
         * \param[in] findex File index (in the energies stored)
         * \param[in] bStoreData boolean indicating whether to store the data
         * \param[in] eTerm  String describing the energy
         * \param[in] eUnit  String describing the energy unit
         */
        EnergyTerm(unsigned int findex,
                   bool         bStoreData,
                   std::string  eTerm,
                   std::string  eUnit);

        //! Return the index in the file to the function type stored here
        unsigned int fileIndex() const { return findex_; }

        //! Return the name corresponding to the energy term
        const std::string name() const { return eTerm_; }

        //! Return the name corresponding to the energy unit
        const std::string unit() const { return eUnit_; }

        /*! \brief
         * Tell the class to store or not to store data
         * \param[in] bStoreData Boolean
         */
        void setStoreData(bool bStoreData) { storeData_ = bStoreData; }

        //! Return the store data variable
        bool storeData() const { return storeData_; }

        //! Is this a true energy or e.g. Temperature
        bool isEner() const { return isEner_; }

        //! Return iterator to begin looping over energy frames
        EnergyFrameIterator begin() const { return ef_.begin(); }

        //! Return iterator to end looping over energy frames
        EnergyFrameIterator end() const { return ef_.end(); }

        /*! \brief
         * Return the energy frame corresponding to a certain step
         * \param[in] iframe The index in the array
         * \return the actual EnergyFrameIterator, or end() if not found
         */
        EnergyFrameIterator findFrame(gmx_int64_t iframe) const;

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
        double average() const { return energy_; }

        //! Return the standard deviation
        double standardDeviation() const { return stddev_; }

        /*! \brief
         * Compute an error estimate based on block averaging.
         * Requires that the energies have been stored.
         * \param[in] nb Number of blocks
         * \param[out] ee the error estimate
         * \return true if an error estimate was computed
         */
        bool errorEstimate(unsigned int nb, real *ee) const;

        /*! \brief
         * Calculate the drift - can only be done when the data is stored.
         * This is done by fitting the data to a line y = ax + b.
         * \param[out] drift The slope of the line (property per time unit)
         * \return true if the result was indeed calculated
         */
        bool drift(real *drift) const;

        //! Return the number of points stored
        gmx_int64_t numFrames() const { return ef_.size(); }

        //! Return the length of the data set in time
        double timeSpan() const { return timeEnd() - timeBegin(); }

        //! Return the begin time of the data set
        double timeBegin() const { return t0_; }

        //! Return the end time of the data set
        double timeEnd() const { return t1_; }

        //! Return the length of the data set in time
        gmx_int64_t numSteps() const { return 1 + (stepEnd() - stepBegin()); }

        //! Return the begin time of the data set
        gmx_int64_t stepBegin() const { return step0_; }

        //! Return the end time of the data set
        gmx_int64_t stepEnd() const { return step1_; }
};

} // namespace energyanalysis

} // namespace gmx
#endif
