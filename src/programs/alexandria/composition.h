/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \defgroup module_alexandria Processing of input files and force fields
 * \ingroup group_preprocessing
 * \brief
 * Provides tools for input processing based on the alexandria force field 
 *
 * The tools in the alexandria directory are under heavy development still.
 * Once finished they will provide processing of small molecules files based
 * on quantum chemistry or other input source through the OpenBabel library.
 * Assigning of atom types, derivation of charges, and generation of 
 * molecular topology files (or objects) is implemented.
 *
 * \author David van der Spoel <david.vanderspoel@gmail.com>
 * \inpublicapi
 * \ingroup module_alexandria
 */
#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <string>
#include <vector>

namespace alexandria
{
/*! \brief
 * enum describing the different composition schemes
 * \inpublicapi
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
    iComp ic_;
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
    CompositionSpec(iComp ic, const char *name, const char *reference, const char *abbreviation) {
        ic_ = ic;
        name_.assign(name); reference_.assign(reference); abbreviation_.assign(abbreviation);
    }
    
    ~CompositionSpec() {};
    
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
    
    CompositionSpecIterator searchCS(iComp ic) { 
        CompositionSpecIterator csi;
        for(csi=beginCS(); (csi<endCS()); ++csi)
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
