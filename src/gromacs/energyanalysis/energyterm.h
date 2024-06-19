/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * \brief
 * Declares gmx::EnergyTerm
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ENERGYTERM_H
#define GMX_ENERGYANALYSIS_ENERGYTERM_H

#include <cstdint>

#include <optional>
#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "energyanalysisframe.h"

namespace gmx
{

//! Typedef for looping over EnergyFrame
using EnergyAnalysisFrameIterator = std::vector<EnergyAnalysisFrame>::const_iterator;

/*! \libinternal
 * \brief
 * Class describing the whole time series of an energy term.
 */
class EnergyTerm
{
private:
    //! Name of the energy term
    std::string energyTerm_;
    //! Unit of this energy
    std::string energyUnit_;
    //! Number of energy terms summed so far
    int64_t numberOfEnergyTerms_ = 0;
    //! First MD step in the analysis
    int64_t firstStep_ = 0;
    //! Last MD step in the analysis added so far
    int64_t lastStep_ = 0;
    //! Index in the energy array in the energy file
    unsigned int indexWithinEnergyFile_;
    //! Best estimate of the average energy so far
    double average_ = 0;
    //! Best estimate of the standard deviation so far
    double standardDeviation_ = 0;
    //! Start time of the analysis
    double startTime_ = 0;
    //! End time of the analysis so far
    double endTime_ = 0;
    //! Boolean whether firstframe has been read
    bool firstFrameRead_ = false;
    //! Boolean indicating whether we are storing data in the vector below
    bool storeData_;
    //! Array of energy frames
    std::vector<EnergyAnalysisFrame> energyAnalysisFrames_;
    //! Total sum of energy
    double totalSumOfEnergy_ = 0;
    //! Total variance of energy
    double totalVarianceOfEnergy_ = 0;
    //! Is the present energy term really an energy?
    bool termIsEnergy_ = false;

public:
    /*! \brief
     * Constructor
     * \param[in] indexWithinEnergyFile File index (in the energies stored)
     * \param[in] bStoreData boolean indicating whether to store the data
     * \param[in] energyTerm  String describing the energy
     * \param[in] energyUnit  String describing the energy unit
     */
    EnergyTerm(unsigned int       indexWithinEnergyFile,
               bool               bStoreData,
               const std::string& energyTerm,
               const std::string& energyUnit);

    //! Return the index in the file to the function type stored here
    unsigned int fileIndex() const { return indexWithinEnergyFile_; }

    //! Return the name corresponding to the energy term
    std::string name() const { return energyTerm_; }

    //! Return the name corresponding to the energy unit
    std::string unit() const { return energyUnit_; }

    /*! \brief
     * Tell the class to store or not to store data
     * \param[in] bStoreData Boolean
     */
    void setStoreData(bool bStoreData) { storeData_ = bStoreData; }

    //! Returns whether any data has been stored from the analyzed frame
    bool storeData() const { return storeData_; }

    //! Is this a true energy or e.g. Temperature
    bool termIsEnergy() const { return termIsEnergy_; }

    //! Return iterator to begin looping over energy frames
    EnergyAnalysisFrameIterator begin() const { return energyAnalysisFrames_.begin(); }

    //! Return iterator to end looping over energy frames
    EnergyAnalysisFrameIterator end() const { return energyAnalysisFrames_.end(); }

    /*! \brief
     * Return the stored energy frame corresponding to a certain input frame.
     * \param[in] frameIndex The frame number to search the storage for.
     * \return the actual EnergyFrameIterator, or end() if not found
     */
    EnergyAnalysisFrameIterator findFrame(int64_t frameIndex) const;

    /*! \brief
     * Add a data frame to this EnergyTerm
     * \param[in] time The time in the simulation
     * \param[in] step The simulation step
     * \param[in] numIntermediateStepsSum The number of intermediate steps for the sums
     * \param[in] energySumOverNumSteps The sum of energies over the last nsum steps
     * \param[in] energyVarianceOverNumSteps The variance of the energies over the last nsum steps
     * \param[in] energyAtTime The energy at this point in time (trajectory)
     */
    void addFrame(double  time,
                  int64_t step,
                  int     numIntermediateStepsSum,
                  double  energySumOverNumSteps,
                  double  energyVarianceOverNumSteps,
                  double  energyAtTime);

    //! Return the average energy
    double average() const { return average_; }

    //! Return the standard deviation
    double standardDeviation() const { return standardDeviation_; }

    /*! \brief
     * Compute an error estimate based on block averaging.
     * Requires that the energies have been stored.
     * \param[in] numBlocks Number of blocks
     * \return Optional error estimate
     */
    std::optional<real> errorEstimate(unsigned int numBlocks) const;

    /*! \brief
     * Calculate the slope of linear fit of energy
     *
     * This energy drift calculation can only be done when the data is stored.
     * This is done by fitting the data to a line y = ax + b.
     * \return Optional slope of the line when calculated
     */
    std::optional<real> slopeOfLinearFit() const;

    //! Return the number of points stored
    int64_t numFrames() const { return energyAnalysisFrames_.size(); }

    //! Return the length of the data set in time
    double timeSpan() const { return timeEnd() - timeBegin(); }

    //! Return the begin time of the data set
    double timeBegin() const { return startTime_; }

    //! Return the end time of the data set
    double timeEnd() const { return endTime_; }

    //! Return the length of the data set in steps
    int64_t numSteps() const { return stepEnd() - stepBegin(); }

    //! Return the begin step of the data set
    int64_t stepBegin() const { return firstStep_; }

    //! Return the end step of the data set
    int64_t stepEnd() const { return lastStep_; }
};

} // namespace gmx
#endif
