/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <string>
#include <vector>

namespace alexandria
{
/*! \brief
 * enum describing the different composition schemes
 * \ingroup module_alexandria
 */
enum iComp {
    iCalexandria, iCbosque, iCmiller
};

/*! \brief
 * Return the name corresponding to this composition type
 *
 * \param ic The composition type
 * \return the name
 */

class CompositionSpec
{
    private:
        // Element of enum
        iComp       ic_;
        //! Name of the composition
        std::string name_;
        //! Literature reference
        std::string reference_;
        //! Abbreviation
        std::string abbreviation_;
    public:
        /*! \brief
         * Constructor
         *
         * \param[in] ic The enum for this composition
         * \param[in] name Name of the composition
         * \param[in] reference Literature reference
         * \param[in] abreviation Short notation
         */
        CompositionSpec(iComp ic, const char *name, const char *reference, const char *abbreviation)
        {
            ic_ = ic;
            name_.assign(name); reference_.assign(reference); abbreviation_.assign(abbreviation);
        }

        //~CompositionSpec() {};

        //! Return the iComp
        iComp iC() { return ic_; }

        //! Return the name
        const char *name() { return name_.c_str(); }

        //! Return the reference
        const char *reference() { return reference_.c_str(); }

        //! Return the abbreviation
        const char *abbreviation() { return abbreviation_.c_str(); }
};

typedef std::vector<CompositionSpec>::iterator CompositionSpecIterator;

class CompositionSpecs
{
    private:
        //! Array of composition specifications
        std::vector<CompositionSpec> cs_;
    public:
        //! Constructor
        CompositionSpecs();

        //! Destructor
        ~CompositionSpecs() {}

        CompositionSpecIterator beginCS() { return cs_.begin(); }

        CompositionSpecIterator endCS() { return cs_.end(); }

        CompositionSpecIterator searchCS(iComp ic)
        {
            CompositionSpecIterator csi;
            for (csi = beginCS(); (csi < endCS()); ++csi)
            {
                if (csi->iC() == ic)
                {
                    break;
                }
            }
            return csi;
        }
};

}

#endif
