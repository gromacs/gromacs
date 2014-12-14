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
 * Declares gmx::energyanalysis::EnergyHelper
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_HELPER_H
#define GMX_ENERGYANALYSIS_HELPER_H

#include <string>
#include <vector>

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/options/options.h"

#include "energyterm.h"

namespace gmx
{

namespace energyanalysis
{

//! Typedef for looping over EnergyTerm
typedef std::vector<EnergyTerm>::const_iterator ConstEnergyTermIterator;

//! Typedef for looping over EnergyTerm
typedef std::vector<EnergyTerm>::iterator EnergyTermIterator;

/*! \libinternal
 * \brief
 * Class for extracting energy terms.
 *
 * This code can store multiple energy terms from energy files or from
 * trajectory next generation files. A number of utilities are
 * present, such as keeping track of the number of molecules.
 *
 * The class can also compute average and standard deviations for energy
 * terms named by a string (function energyTerm).
 *
 * In order to check convergence
 * of the simulations a plot of the variance can be generated.
 *
 * Raw energies as well as statistics can be printed to a file.
 */
class EnergyTermContainer
{
    private:
        //! Should we store the data in memory?
        bool                    storeData_;

        //! Number of molecules
        int                     nMol_;

        //! The energy data should be stored here when bRead_ is set
        std::vector<EnergyTerm> et_;

        //! Number of blocks for error analysis
        int                     nBlocks_;

    public:
        //! Constructor
        EnergyTermContainer()
            : storeData_(false), nMol_(1), nBlocks_(5)
        {};

        //! Additional generic command line options
        void initOptions(IOptionsContainer *options);

        /*! \brief
         * Get the numbers of blocks for error estimates
         * \return Number
         */
        unsigned int nBlocks() const { return nBlocks_; }

        //! Return the number of molecules
        int nMol() const { return nMol_; }

        /*! \brief
         * Tell the class to store or not to store data
         * \param[in] storeData Boolean
         */
        void setStoreData(bool storeData);

        //! Return the store data variable
        bool storeData() { return storeData_; }

        //! Store all data from one frame
        void addFrame(t_enxframe *fr);

        //! Get start iterator
        EnergyTermIterator begin() { return et_.begin(); }

        //! Get end iterator
        EnergyTermIterator end() { return et_.end(); }

        //! Get start iterator
        ConstEnergyTermIterator begin() const { return et_.begin(); }

        //! Get end iterator
        ConstEnergyTermIterator end() const { return et_.end(); }

        /*! \brief
         * Search a specific iterator
         * \param[in] findex Index in energy file
         * \return the iterator or end()
         */
        EnergyTermIterator etSearch(unsigned int findex);

        /*! \brief
         * Search a specific iterator
         * \param[in] eTerm The name of the energy term
         * \return the iterator or end()
         */
        EnergyTermIterator etSearch(std::string eTerm);

        //! Add energy term
        void addEnergyTerm(EnergyTerm et) { et_.push_back(et); }

        //! Return number of energy terms
        unsigned int nEnergyTerm() { return et_.size(); }

        /*! \brief
         * Extract a single energy term and its standard deviation.
         *
         * \param[in] ftype GROMACS internal function type
         * \param[out] e  The average energy
         * \param[out] stddev The standard deviation
         * \return true if succesfull, false otherwise
         */
        bool energyTerm(unsigned int ftype, double *e, double *stddev);

        /*! \brief
         * Extract a single energy term and its standard deviation.
         *
         * \param[in] term Name of the energy term
         * \param[out] e  The average energy
         * \param[out] stddev The standard deviation
         * \return true if succesfull, false otherwise
         */
        bool energyTerm(std::string term, double *e, double *stddev);
};


/*! \brief
 * Print the average and other statistics for all energy terms.
 * \param[in] fp      File pointer to print to
 * \param[in] eBegin  Energy term begin iteator
 * \param[in] eEnd    Energy term end iteator
 * \param[in] nBlocks Number of blocks to divide trajectory into
 */
void printStatistics(FILE                   *fp,
                     ConstEnergyTermIterator eBegin,
                     ConstEnergyTermIterator eEnd,
                     unsigned int            nBlocks);

/*! \brief
 * Print the legend for all the stored energy terms to a file pointer.
 * \param[in] fp      The file pointer opened previously.
 * \param[in] eBegin  Energy term begin iteator
 * \param[in] eEnd    Energy term end iteator
 * \param[in] oenv    GROMACS output environment
 */
void printXvgLegend(FILE                   *fp,
                    ConstEnergyTermIterator eBegin,
                    ConstEnergyTermIterator eEnd,
                    const OutputEnvPointer &oenv);

/*! \brief
 * Print the energies to an xvg file. Can only be done if the energies
 * have been stored, by calling setStoreData(true) prior to reading the
 * file.
 * \param[in] outputFile file name for the xvg file to be written.
 * \param[in] eBegin     Energy term begin iteator
 * \param[in] eEnd       Energy term end iteator
 * \param[in] bDouble    Whether or not to print in double precision
 * \param[in] oenv       GROMACS output environment
 */
void printEnergies(const std::string      &outputFile,
                   ConstEnergyTermIterator eBegin,
                   ConstEnergyTermIterator eEnd,
                   bool                    bDouble,
                   const OutputEnvPointer &oenv);

/*! \brief
 * Construct an y-axis with units as GROMACS likes it
 * \param[in]  eBegin Energy term begin iteator
 * \param[in]  eEnd   Energy term end iteator
 * \param[out] yaxis  The string with units
 */
void yAxis(ConstEnergyTermIterator  eBegin,
           ConstEnergyTermIterator  eEnd,
           std::string             *yaxis);

} // namespace energyanalysis

} // namespace gmx

#endif
