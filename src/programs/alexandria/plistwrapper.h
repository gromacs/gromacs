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
#ifndef PLISTWRAPPER_H
#define PLISTWRAPPER_H

//#include <algorithm>
#include <vector>

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"

namespace alexandria
{
//! Interaction type
enum InteractionType
{
    eitBONDS              = 0,
    eitANGLES             = 1,
    eitLINEAR_ANGLES      = 2,
    eitPROPER_DIHEDRALS   = 3,
    eitIMPROPER_DIHEDRALS = 4,
    eitVDW                = 5,
    eitLJ14               = 6,
    eitPOLARIZATION       = 7,
    eitCONSTR             = 8,
    eitVSITE2             = 9,
    eitNR                 = 10
};

using ParamIterator = typename std::vector<t_param>::iterator;

//! Cleaner version of plist array
class PlistWrapper
{
    private:
        //! Function type
        int                  ftype_;
        //! Interaction type
        InteractionType      itype_;
        //! Array of parameters
        std::vector<t_param> p_;
    public:
        //! Constructor
        PlistWrapper(InteractionType itype,
                     int             ftype) : ftype_(ftype), itype_(itype) {}

        //! Add one parameter
        void addParam(t_param p) { p_.push_back(p); }

        //! Return the function type
        int getFtype() const { return ftype_; }

        //! Return the interaction type
        InteractionType getItype() const { return itype_; }

        //! Update the function type
        void setFtype(int ftype) { ftype_ = ftype; }

        //! Return the interaction type
        InteractionType interactionType() const { return itype_; }

        //! Loop over parameters
        ParamIterator beginParam() { return p_.begin(); }

        //! Loop over parameters
        ParamIterator endParam() { return p_.end(); }

        //! Remove one parameter from the array and return array for next
        ParamIterator eraseParam(ParamIterator p) { return p_.erase(p); }

        //! Remove all parameters
        void eraseParams() { p_.clear(); }

        //! Return number of parameters
        unsigned int nParam() { return p_.size(); }
};

//! Another utility typedef for a looper
using  PlistWrapperIterator = typename std::vector<PlistWrapper>::iterator;

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype);

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, InteractionType itype);

unsigned int CountPlist(std::vector<PlistWrapper> &plist, int ftype);

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int                  ftype,
                        const t_param             &p);

void delete_params(std::vector<PlistWrapper> &plist_,
                   const int                  ftype,
                   const int                  alist[]);

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        int                        ftype,
                        InteractionType            itype,
                        const t_param             &p);
}

#endif
