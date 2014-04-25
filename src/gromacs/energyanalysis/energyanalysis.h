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
#ifndef GMX_ENERGYANALYSIS_H
#define GMX_ENERGYANALYSIS_H

#include <string>
#include <vector>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/fileio/enxio.h"
#include "energyterm.h"

namespace gmx
{
/*! \brief
 * Class for extracting energy terms from energy files or from
 * trajectory next generation files. In order to check convergence
 * of the simulations a plot of the variance can be generated.
 */
class EnergyAnalysis
{
    private:
        //! Currently/most recently used energy filename
        std::string eneFile_;

        //! Has the file been read?
        bool bRead_;

        //! Should we write in double precision?
        bool bDouble_;

        //! Should we store the data in memory?
        bool bStoreData_;

        //! Number of molecules
        unsigned int          nMol_;

        //! Output file name
        std::string outputFile_;

        //! Output data
        output_env_t oenv_;

        //! The energy data should be stored here when bRead_ is set
        std::vector<EnergyTerm> et_;

        //! Number of blocks for error estimate
        unsigned int nB_;

    public:
        //! Default constructor
        EnergyAnalysis()
        {
            bRead_      = false;
            bDouble_    = false;
            bStoreData_ = false;
            oenv_       = NULL;
            nMol_       = 1;
            nB_         = 1;
        };

        //! Destructor
        ~EnergyAnalysis() {};

        //! Set the energy file name
        void setEneFile(std::string eneFile) { eneFile_ = eneFile; }

        //! Get the energy file name
        std::string getEneFile() { return eneFile_; }

        //! Set the output file name
        void setOutputFile(std::string outputFile) { outputFile_ = outputFile; }

        //! Get the output file name
        std::string getOutputFile() { return outputFile_; }

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
        unsigned int getNblocks() { return nB_; }

        //! \return the number of molecules
        int getNmol() { return nMol_; }

        /*! \brief
         * Does the initiation of the analysis of the file
         * \param[in] nre Number of energy terms in the file
         * \param[in] enm Names of the energy terms etc.
         * \return true if OK, false otherwise
         */
        virtual bool initAnalysis(int nre, gmx_enxnm_t enm[]) = 0;

        /*! \brief
         * Analyse one frame and stores the results in memory
         * \param[in] fr The energy data frame
         * \return true if OK, false otherwise
         */
        virtual bool addAnalysisFrame(t_enxframe *fr) = 0;

        //! Finalize reading
        virtual bool finalizeAnalysis() = 0;

        //! Do the reading by calling the above three routines
        bool readEneFile();

        /*! \brief
         * Tell the class to store or not to store data
         * \param[in] bStoreData Boolean
         */
        void setStoreData(bool bStoreData);

        //! Return the store data variable
        bool getStoreData() { return bStoreData_; }

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
        output_env_t getOutputEnvironment() { return oenv_; }

        //! Return number of energy terms (will be 0 before readEneFile() is called)
        unsigned int nEnergyTerm() { return et_.size(); }

        //! Set the double precision flag
        void setDoublePrecision(bool bDouble) { bDouble_ = bDouble; }

        //! Get the double precision flag
        bool getDoublePrecision() { return bDouble_; }

        /*! \brief
         * Extract a single energy term. Only when this routine is called is
         * the file read and, if requested, the fluctuation convergence is stored.
         * \param[in] ftype GROMACS internal function type
         * \param[out] e  The average energy
         * \param[out] stddev The standard deviation
         * \return true if succesfull, false otherwise
         */
        bool getEnergyTerm(unsigned int ftype, double *e, double *stddev);

        /*! \brief
         * Extract a single energy term. Only when this routine is called is
         * the file read and, if requested, the fluctuation convergence is stored.
         * \param[in] term Name of the energy term
         * \param[out] e  The average energy
         * \param[out] stddev The standard deviation
         * \return true if succesfull, false otherwise
         */
        bool getEnergyTerm(const char *term, double *e, double *stddev);

        //! Sum all the energy terms and delete the original data sets
        void sumEnergies();

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
        void printEnergies();

        //! View the output file(s)
        void viewOutput();
};

} // namespace gmx

#endif
