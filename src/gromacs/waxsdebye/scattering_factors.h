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
#ifndef _scattering_factors_h
#define _scattering_factors_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <vector>
#include <libxml/tree.h>
#include <libxml/parser.h>
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

/*! \brief
 * Base class for storing scattering factor parameters for an atom type.
 */
class ScatteringFactor
{
    protected:
        //! Residue name
        std::string residue_;
        //! Atom name
        std::string atom_;
        //! Type corresponding to atom
        int         type_;
    public:
        //! Constructor
        ScatteringFactor() { type_ = 0; }

        //! Destructor
        virtual ~ScatteringFactor() {};

        //! Returns the type of this atom
        int type() const { return type_; }

        //! Returns the residue  as a string
        const std::string &residue() const { return residue_; }

        //! Returns the atom symbol as a string
        const std::string &atom() const { return atom_; }

        /*! \brief
         * Virtual function to compute the scattering for a given q
         * \param[in] q scattering vector
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double q) const = 0;

        /*! \brief
         * Virtual function to compute the scattering for a given q
         * \param[in] theta scattering angle
         * \param[in] lambda wave length, should be > 0
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double theta, double lambda) const = 0;

        /*! \brief
         * Function to load data from an array retrieved from an xml structure
         * \param[in] tree the libxml2 data structure
         */
        virtual void loadScatteringFactor(xmlNodePtr tree) = 0;

        /*! \brief
         * Function to store data in xml structure
         * \param[in] parent xml data structure
         */
        virtual void storeScatteringFactor(xmlNodePtr parent) = 0;
};

//! Typedef for encapsulating scattering factor array
typedef gmx_unique_ptr<ScatteringFactor>::type ScatteringFactorPointer;

//! Enumerated type to distinguish the parameters
enum esfType {
    esfCromerMann, esfFourier
};

/*! \brief
 * Class for loading and storing scattering factor parameters for several
 * atom types, as well as calculating f(q) for an atom. Note that the
 * elements do not have to chemical elements, they may be e.g. coarse
 * grained particles or entire amino acids.
 */
class ScatteringFactorTable
{
    private:
        // Force field name
        std::string ForceField_;
        //! Reference to papers
        std::string Reference_;
        // Type of the elements (Fourier or CromerMann)
        esfType     esftype_;
        // Is the displaced solvent included in the parameters?
        bool        bDisplacedSolvent_;
        //! Vector of length number of elements in table.
        std::vector<ScatteringFactorPointer> cm_;
    public:
        //! Constructor
        ScatteringFactorTable();

        //! Destructor
        ~ScatteringFactorTable() {}

        //! Begin iterator
        std::vector<ScatteringFactorPointer>::iterator beginSFP() { return cm_.begin(); }

        //! End iterator
        std::vector<ScatteringFactorPointer>::iterator endSFP() { return cm_.end(); }

        //! Add table entry
        void addSFP(ScatteringFactorPointer sfp) { cm_.push_back(gmx::move(sfp)); }
        /*! \brief
         * Read table from a file
         *
         * \param[in]  datafile  File containing the scattering factors
         * \return true if reading went fine, false otherwise
         */
        bool read(const char *datafile);

        /*! \brief
         * Write table to a file
         *
         * \param[in] datafile  File to write the scattering factors to
         * \return true if writing went fine, false otherwise
         */
        bool write(const char *datafile);

        //! Set the parameter type
        void setSfType(esfType esftype) { esftype_ = esftype; }

        //! Set the force field name
        void setForceField(const char *force_field) { ForceField_.assign(force_field); }

        //! Set the references
        void setReference(const char *reference) { Reference_.assign(reference); }

        //! Set the solvent layer
        void setDisplacedSolvent(bool bDisplacedSolvent) { bDisplacedSolvent_ = bDisplacedSolvent; }

        //! Return the parameter type
        esfType getSfType() { return esftype_; }

        //! Return the force field name
        std::string getForceField() { return ForceField_; }

        //! Return the references
        std::string getReference() { return Reference_; }

        //! Return whether a solvent layer was used
        bool getDisplacedSolvent() { return bDisplacedSolvent_; }

        //! Return length of the table
        unsigned int size() const { return cm_.size(); }

        /*! \brief
         * \param i index in table
         * \return the type corresponding to index
         */
        int type(int i) const { return cm_[i]->type(); }

        /*! \brief
         * Compute scattering factor
         * \param[in] type  Number of the element looked for
         * \param[in] q     Scattering vector
         * \return scattering factor
         */
        double computeScatteringFactor(int type, double q) const;

        /*! \brief
         * Compute scattering factor
         * \param[in] residue  Residue name
         * \param[in] atom     Atom name
         * \param[in] q        Scattering vector
         * \return scattering factor
         */
        double computeScatteringFactor(const std::string &residue,
                                       const std::string &atom,
                                       double             q) const;

        //! Return largest type in the table
        int maxType() const;
};

}


#endif /* _scattering_factors check */
