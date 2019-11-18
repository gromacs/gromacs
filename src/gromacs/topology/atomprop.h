/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_ATOMPROP_H
#define GMX_TOPOLOGY_ATOMPROP_H

#include <string>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

enum
{
    epropMass,
    epropVDW,
    epropDGsol,
    epropElectroneg,
    epropElement,
    epropNR
};

struct AtomProperty;
class ResidueType;
/*! \brief
 * Holds all the atom property information loaded.
 */
class AtomProperties
{
public:
    //! Default constructor.
    AtomProperties();
    //! Default destructor
    ~AtomProperties();

    /*! \brief
     * Get element string from atom number.
     *
     * \param[in] atomNumber Atomnumber to check.
     * \returns Name of the element.
     *
     * \todo This should be made const once the lazy
     * implementation is done properly for the class.
     */
    std::string elementFromAtomNumber(int atomNumber);
    /*! \brief
     * Get atom number from element string.
     *
     * \param[in] element Name of element.
     * \returns AtomNumber that was being looked for.
     *
     * \todo This should be made const once the lazy
     * implementation is done properly for the class.
     */
    int atomNumberFromElement(const char* element);
    /*! \brief
     * Set atom property based on atomname.
     *
     * Extract a \p value from the database. Returns true
     * if this is successful, or false if not. Sets default value
     * in the later case. The first time this function is called
     * for this property the database will be initialized.
     *
     * \param[in] eprop Property to set.
     * \param[in] residueName Residue name for entry.
     * \param[in] atomName Atom name for entry.
     * \param[out] value New value to set or default.
     * \returns If the operation has been succesful.
     */
    bool setAtomProperty(int eprop, const std::string& residueName, const std::string& atomName, real* value);
    /*! \brief
     * Get handle to property.
     *
     * \param[in] eprop Which property we need a handle to.
     * \returns Pointer to property entry.
     */
    AtomProperty* prop(int eprop);
    //! Get handle to residuetype library.
    ResidueType* restype();

private:
    //! Implementation pointer.
    class Impl;

    gmx::PrivateImplPointer<Impl> impl_;
};

#endif
