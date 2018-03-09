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
 
 
#ifndef MOLSELECT_H
#define MOLSELECT_H

#include <algorithm>
#include <string>
#include <vector>

#include "molprop.h"
#include "poldata.h"

enum iMolSelect {
    imsTrain, imsTest, imsIgnore, imsUnknown, imsNR
};

const char *iMolSelectName(iMolSelect ims);

namespace alexandria
{

class IMolSelect 
{
    private:
        std::string iupac_;
        iMolSelect  status_;
        int         index_;
        
    public:
        IMolSelect(const std::string &iupac, iMolSelect status, int index)
            :
                iupac_(iupac), 
                status_(status), 
                index_(index) 
            {}

        const std::string &iupac() const { return iupac_; }
        
        iMolSelect status() const { return status_; }
        
        int index() const { return index_; }
};

class MolSelect
{
    private:
    
        std::vector<IMolSelect> ims_;

    public:
    
        MolSelect() {};

        void read(const char *filename);
        
        size_t nMol() const { return ims_.size(); }
        
        iMolSelect status(const std::string &iupac) const;
        
        int index(const std::string &iupac) const;
        
        int count(iMolSelect ims) const
        {
            return std::count_if(ims_.begin(), ims_.end(),
                                 [ims](IMolSelect const &i)
                                 { return i.status() == ims; });
        }
};

} // namespace

#endif
