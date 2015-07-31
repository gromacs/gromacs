/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
#include "gromacs/options/options.h"

#include "energyterm.h"

namespace gmx
{

namespace energyanalysis
{

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
class EnergyHelper
{
    private:
        //! Should we write in double precision?
        bool bDouble_;

        //! Should we store the data in memory?
        bool bStoreData_;

        //! Number of molecules
        int          nMol_;

        //! Output data
        output_env_t oenv_;

        //! The energy data should be stored here when bRead_ is set
        std::vector<EnergyTerm> et_;

        //! Number of blocks for error estimate
        unsigned int nB_;

    public:
        //! Constructor
        EnergyHelper()
        {
            bDouble_    = false;
            bStoreData_ = false;
            oenv_       = NULL;
            nMol_       = 1;
            nB_         = 1;
        };

        //! Additional generic command line options
        void initOptions(IOptionsContainer *options);

        //! Set the number of molecules to nmol
        void setNmol(int nMol) { nMol_ = nMol; }

        /*! \brief
         * Set the numbers of blocks for error estimates
         * \param[in] nb Number
         */
        void setNblocks(unsigned int nb) { nB_ = nb; }

        /*! \brief
         * Get the numbers of blocks for error estimates
         * \return Number
         */
        unsigned int nBlocks() { return nB_; }

        //! Return the number of molecules
        int nMol() { return nMol_; }

        /*! \brief
         * Tell the class to store or not to store data
         * \param[in] bStoreData Boolean
         */
        void setStoreData(bool bStoreData);

        //! Return the store data variable
        bool storeData() { return bStoreData_; }

        //! Get start iterator
        EnergyTermIterator etBegin() { return et_.begin(); }

        //! Get end iterator
        EnergyTermIterator etEnd() { return et_.end(); }

        /*! \brief
         * Search specific iterator
         * \param[in] findex Index in energy file
         * \return the iterator
         */
        EnergyTermIterator etSearch(unsigned int findex);

        /*! \brief
         * Search specific iterator
         * \param[in] eTerm The name of the energy term
         * \return the iterator
         */
        EnergyTermIterator etSearch(std::string eTerm);

        //! Add energy term
        void addEnergyTerm(EnergyTerm et) { et_.push_back(et); }

        //! Set output environment
        void setOutputEnvironment(output_env_t oenv) { oenv_ = oenv; }

        //! Get output environment
        output_env_t outputEnvironment() { return oenv_; }

        //! Return number of energy terms
        unsigned int nEnergyTerm() { return et_.size(); }

        //! Set the double precision flag
        void setDoublePrecision(bool bDouble) { bDouble_ = bDouble; }

        //! Get the double precision flag
        bool doublePrecision() { return bDouble_; }

        /*! \brief
         * Extract a single energy term.
         *
         * Only when this routine is called is
         * the file read and, if requested, the fluctuation convergence is stored.
         * \param[in] ftype GROMACS internal function type
         * \param[out] e  The average energy
         * \param[out] stddev The standard deviation
         * \return true if succesfull, false otherwise
         */
        bool energyTerm(unsigned int ftype, double *e, double *stddev);

        /*! \brief
         * Extract a single energy term. Only when this routine is called is
         * the file read and, if requested, the fluctuation convergence is stored.
         * \param[in] term Name of the energy term
         * \param[out] e  The average energy
         * \param[out] stddev The standard deviation
         * \return true if succesfull, false otherwise
         */
        bool energyTerm(std::string term, double *e, double *stddev);

        /*! \brief
         * Print the average and other statistics for all energy terms.
         * \param[in] fp File pointer to print to
         */
        void printStatistics(FILE *fp);

        /*! \brief
         * Print the legend for all the stored energy terms to a file pointer.
         * \param[in] fp The file pointer opened previously.
         */
        void printXvgLegend(FILE *fp);

        /*! \brief
         * Print the energies to an xvg file. Can only be done if the energies
         * have been stored, by calling setStoreData(true) prior to reading the
         * file.
         */
        void printEnergies(std::string outputFile);
};

}

} // namespace gmx

#endif
