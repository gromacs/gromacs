/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GAUSS_IO_H
#define GAUSS_IO_H

#include "molprop.h"

namespace alexandria
{

/*! \brief
 * Helper class for extracting atomization energies from Gaussian calcs.
 *
 * This contains data in order to extract enthalpy of
 * formation and Gibbs energy of formation from Thermochemistry methods
 * in Gaussian.
 *
 * \inlibraryapi
 * \ingroup module_alexandria
 */
class GaussAtomPropVal
{
    private:
        std::string _element, _method, _desc;
        real        _temp;
        real        _value;
    public:
        //! Constructor initiating all the values stored
        GaussAtomPropVal(std::string element, std::string method, std::string desc,
                         real temp, real value)
        {
            _element = element; _method = method; _desc = desc;
            _temp    = temp; _value = value;
        }

        //! Default destructor
        ~GaussAtomPropVal() {}

        //! Return element name
        std::string getElement() { return _element; }

        //! Return the method used
        std::string getMethod() { return _method; }

        //! Return a description of the type of value stored
        std::string getDesc() { return _desc; }

        //! Return the temperature
        real getTemp() { return _temp; }

        //! Return the actual value
        real getValue() { return _value; }
};

/*! \brief
 * Class for extracting atomization energies from Gaussian calcs.
 *
 * This contains data in order to extract enthalpy of
 * formation and Gibbs energy of formation from Thermochemistry methods
 * in Gaussian. The method and description must match that in a library file.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class GaussAtomProp
{
    private:
        std::vector<GaussAtomPropVal> _gapv;
    public:
        //! Default constructor
        GaussAtomProp();

        //! Default destructor
        ~GaussAtomProp() {}

        /*! \brief
         * Look up the value corresponding to input variables
         *
         * \param[in] element  From the periodic table
         * \param[in] method   Method used, e.g. G3, exp
         * \param[in] desc     Which type of name is
         * \param[in] temp     Temperature
         * \param[out] value   The energy value
         * \return 1 on success, 0 otherwise
         * \ingroup module_alexandria
         */
        int getValue(const char *element,
                     const char *method,
                     const char *desc,
                     double      temp,
                     double     *value);
};
}

/*! \brief
 * Read a Gaussian log file either using home grown methods or using OpenBabel
 *
 *
 * \param[in] g98        The gaussian log file, or in case OpenBabel is used anything
 *                       that can be read by OpenBabel
 * \param[out] mpt       The MolProp
 * \param[in] molnm      Molecule name to override the one from the filename [ maybe NULL ]
 * \param[in] iupac      IUPAC name to override the one from the filename [ maybe NULL ]
 * \param[in] conformation  Conformation the molecule is in [ maybe NULL ]
 * \param[in] basisset   Basis set used for the calculation [ maybe NULL ]
 * \param[in] maxpot     Maximum number of electrostatic potential data points to store
 * \param[in] nsymm      Symmetry number for this molecule. If zero it will be detected from
 *                       the input.
 * \param[in] forcefield One of the force fields supported by OpenBabel used for atomtypes
 * \ingroup module_alexandria
 */
void ReadGauss(const char          *g98,
               alexandria::MolProp &mp,
               const char          *molnm,
               const char          *iupac,
               const char          *conf,
               const char          *basis,
               int                  maxpot,
               int                  nsymm,
               const char          *forcefield,
               const char          *jobtype);

#endif
