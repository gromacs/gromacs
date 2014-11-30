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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef PLISTWRAPPER_H
#define PLISTWRAPPER_H

//#include <algorithm>
#include <vector>

#include "gromacs/gmxpreprocess/grompp-impl.h"

namespace alexandria
{
//! Utility typedef
typedef std::vector<t_param>::iterator ParamIterator;
//! Cleaner version of plist array
class PlistWrapper {
 private:
    int                  ftype_;
    std::vector<t_param> p_;
 public:
    PlistWrapper(int ftype) { setFtype(ftype); }
    ~PlistWrapper() {};
    void addParam(t_param p) { p_.push_back(p); }
    int getFtype() { return ftype_; }
    void setFtype(int ftype) { ftype_ = ftype; }
    ParamIterator beginParam() { return p_.begin(); }
    ParamIterator endParam() { return p_.end(); }
    ParamIterator eraseParam(ParamIterator p) { return p_.erase(p); }
    ParamIterator eraseParams() { return p_.erase(p_.begin(), p_.end()); }
    unsigned int nParam() { return p_.size(); }
};
    
//! Another utility
typedef std::vector<PlistWrapper>::iterator PlistWrapperIterator;
 
PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype);

unsigned int CountPlist(std::vector<PlistWrapper> &plist, int ftype); 
 
void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int ftype,
                        const t_param &p);
                        
void delete_params(std::vector<PlistWrapper> &plist_,
                   const int ftype, 
                   const int alist[]);
                   
void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int ftype,
                        const t_param &p);
}

#endif
