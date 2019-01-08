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

#include <memory>
#include <string>
#include <vector>

#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

enum {
    epropMass, epropVDW, epropDGsol, epropElectroneg, epropElement,
    epropNR
};

//! Basic entries in atomproperties.
struct BaseEntry {
    //! Default constructor.
    BaseEntry(const std::string &aName, const std::string &rName)
        : atomName(aName), residueName(rName), isAvailable(false), value(0.0)
    {}
    //! Name for atom.
    std::string atomName;
    //! Name for residue.
    std::string residueName;
    //! Is property available.
    bool        isAvailable;
    //! Value set for property.
    real        value;
};

//! Conglomeration of basic entries.
struct AtomProperty {
    //! Has property been set.
    bool                   isSet = false;
    //! Database the property is coming from.
    std::string            db;
    //! Default value for property.
    double                 def = 0.0;
    //! Basic entries for properties.
    std::vector<BaseEntry> entry;
};
//! Datastructure containing all atom properties.
struct AtomProperties {
    //! Default constructor.
    AtomProperties() : bWarned(false), bWarnVDW(false), restype(nullptr)
    {
        gmx_residuetype_init(&restype);
    }
    //! Need destructor to clean up residuetype.
    ~AtomProperties()
    {
        gmx_residuetype_destroy(restype);
    }
    //! Has user been warned about error.
    bool               bWarned;
    //! Has user been warned about vdW error.
    bool               bWarnVDW;
    //! The different atom properties.
    AtomProperty       prop[epropNR];
    //! The residue types.
    gmx_residuetype_t *restype;
};


//! Convenvience definition for pointer to AtomProperties.
using AtomPropertiesPtr = std::unique_ptr<AtomProperties>;

/*! \brief
 * Get name of the element for \p atomnumber.
 *
 * \param[in] aps Atom properties data.
 * \param[in] atomNumber Atomnumber to check.
 * \returns Name of the element.
 */
std::string elementFromAtomNumber(AtomProperties *aps, int atomNumber);

/*! \brief
 * Get atomnumber from \p element name.
 *
 * \param[in] aps Atom properties data.
 * \param[in] element Name of element.
 * \returns AtomNumber that was being looked for.
 */
int atomNumberFromElement(AtomProperties *aps, const char *element);

/*! \brief
 * Set atom property based on atomname.
 *
 * \param[in] aps Atom properties data.
 * \param[in] eprop Property to set.
 * \param[in] residueName Residue name for entry.
 * \param[in] atomName Atom name for entry.
 * \param[in,out] value New value to set or default.
 * \returns If the operation has been succesful.
 */
bool setAtomProperty(AtomProperties    *aps,
                     int                eprop,
                     const std::string &residueName,
                     const std::string &atomName,
                     real              *value);
#endif
