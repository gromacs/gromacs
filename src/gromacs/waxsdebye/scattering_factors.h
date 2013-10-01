/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

/*! \brief
 * Base class for storing scattering factor parameters for an atom type.
 */
class ScatteringFactor
{
    protected:
        //! Element name
        std::string element_;
        //! Atomic number corresponding to element
        int         atomic_number_;
    public:
        //! Constructor
        ScatteringFactor() { atomic_number_ = 0; }

        //! Destructor
        virtual ~ScatteringFactor() {};

        //! Returns the atomic number of this element
        int atomic_number() { return atomic_number_; }

        //! Returns the element symbol as a string
        std::string const element() { return element_; }

        /*! \brief
         * Virtual function to compute the scattering for a given q
         * \param[in] q scattering vector
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double q) = 0;

        /*! \brief
         * Virtual function to compute the scattering for a given q
         * \param[in] theta scattering angle
         * \param[in] lambda wave length, should be > 0
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double theta, double lambda) = 0;
};

//! Typedef for encapsulating scattering factor array
typedef gmx_unique_ptr<ScatteringFactor>::type ScatteringFactorPointer;

/*! \brief
 * Class for loading and storing scattering factor parameters for several atom types,
 * as well as calculating f(q) for an atom. Note that the elements do not have to
 * chemical elements, they may be e.g. coarse grained particles or entire amino acids.
 */
class ScatteringFactorTable
{
    private:
        // Vectors of length number of elements in table
        std::vector<ScatteringFactorPointer> cm_;
    public:
        /*! \brief
         * Constructor
         * \param[in] datafile  File containing the scattering factors
         */
        explicit ScatteringFactorTable(const char *datafile);

        //! Destructor
        ~ScatteringFactorTable() {}

        //! Return length of the table
        unsigned int size() const { return cm_.size(); }

        /*! \brief
         * \param i index in table
         * \return the atomic number corresponding to index
         */
        int atomic_number(int i) const { return cm_[i]->atomic_number(); }

        /*! \brief
         * Compute scattering factor
         * \param[in] atomic_number  Atomic Number of the element looked for
         * \param[in] q              Scattering vector
         * \return scattering factor
         */
        double computeScatteringFactor(int atomic_number, double q);

        /*! \brief
         * Compute scattering factor
         * \param[in] element  Element name of the element looked for
         * \param[in] q        Scattering vector
         * \return scattering factor
         */
        double computeScatteringFactor(const std::string &element, double q);

        /*! \brief
         * Compute scattering factor
         * \param[in] element  Element name of the element looked for
         * \param[in] q        Scattering vector
         * \return scattering factor
         */
        double computeScatteringFactor(const char *element, double q)
        {
            std::string s(element);
            return computeScatteringFactor(s, q);
        }

        //! Return largest atomic number in the table
        int max_atomic_number();
};
}

#endif /* _scattering_factors check */
