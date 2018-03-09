/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
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
    eitVSITE3FAD          = 10,
    eitVSITE3OUT          = 11,
    eitNR                 = 12
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
