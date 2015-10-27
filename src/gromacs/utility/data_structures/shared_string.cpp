/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \file
    \ingroup module_utility
    \internal

   \brief
    implementation file for class shared_string

   \author
    R. Thomas Ullmann <tullman@gwdg.de>

   \copyright
    GROMACS license

   \date Mar 2015
 */

#include "gmxpre.h"

#include "shared_string.h"

namespace gmx
{

// replace content of this shared_string and all other instances referring to the same string with the content of s
void shared_string::replaceAll(const std::string &s)
{
    if (strPtr_ == nullptr)
    {
        strPtr_ = getSymtab().nonconstRenameEntry("", s.c_str());
    }
    else
    {
        char * tmpp = (*strPtr_)->getPtr();
        if (tmpp == nullptr)
        {
            strPtr_ = getSymtab().nonconstRenameEntry("", s.c_str());
        }
        else
        {
            // reassigning the shared_ptr can save memory if *this is the
            // only shared_string instance referencing the current string
            strPtr_ = getSymtab().nonconstRenameEntry(tmpp, s.c_str());
        }
    }
}

// append contents of gmx::shared_string rhs to shared_string
shared_string &shared_string::operator+=(const shared_string &rhs)
{
    std::string tmpstr(std::string((*strPtr_)->getPtr()) + std::string((*rhs.strPtr_)->getPtr()));
    strPtr_.reset();
    strPtr_ = getSymtab().nonconstInsertEntry(tmpstr);
    return *this;
}

// append contents of std::string rhs to shared_string
shared_string &shared_string::operator+=(const std::string &rhs)
{
    // by not directly calling insert_entry, it is ensured that the old entry position
    // is already free for the new version if *this was the only shared_ptr instance
    // referencing the old entry
    std::string tmpstr(std::string((*strPtr_)->getPtr()) + rhs);
    strPtr_.reset();
    strPtr_ = getSymtab().nonconstInsertEntry(tmpstr);
    return *this;
}

// defragment the symbol table, remove unused entries and unused buffers
void shared_string::clean()
{
    return getSymtab().clean();
}

// assign left-hand-side shared_string individually to right-hand-side string
shared_string &shared_string::operator=(const std::string &rhs)
{
    strPtr_ = getSymtab().nonconstInsertEntry(rhs);
    return *this;
}
// assign left-hand-side shared_string individually to right-hand-side shared string
shared_string &shared_string::operator=(const char * const rhs)
{
    strPtr_ = getSymtab().nonconstInsertEntry(rhs);
    return *this;
}

// read a string from an istream into the shared_string -- shared_string << istr
std::istream &shared_string::operator<<(std::istream &ist)
{
    return read(ist);
}

// free the symbol table, leaving the entries themselves intact
void shared_string::freeSymtab()
{
    return getSymtab().freeSymtab(false);
}

// print statistics on the symbol table
void shared_string::getStringStats(size_t *ninst, size_t *nptr)
{
    getSymtab().getStringStats(std_str().c_str(), ninst, nptr);
}

// get the amount of memory currently occupied by the symbol table (in bytes)
size_t shared_string::getTotalMem()
{
    return getSymtab().getTotalMem();
}

// get the number of strings currently stored in the symbol table
size_t shared_string::getNr()
{
    return getSymtab().getNr();
}

// get the number of buffers currently used by the symbol table
size_t shared_string::getNbuf()
{
    return getSymtab().getNbuf();
}

// print debug information regarding this shared_string and other shared_string instances referencing the same string content
void shared_string::debugStats()
{
    if (strPtr_)
    {
        printf("  1st level is OK.\n");
        if (*strPtr_)
        {
            printf("    2nd level is OK. use count %li", strPtr_.use_count());
            printf("    (*strPtr_)->getRef()       : %li\n", (*strPtr_)->getRef());
            printf("    (*strPtr_)->getPtr()       : %p\n", static_cast<void*>((*strPtr_)->getPtr()));
            if ((*strPtr_)->getPtr())
            {
                printf("    (*strPtr_)->getPtr() std::string: \"%s\"\n", (*strPtr_)->getPtr());
            }
        }
        else
        {
            printf("    2nd level is NOT OK.\n");
        }
    }
    else
    {
        printf("  1st level is NOT OK.\n");
    }
    size_t ninst = 0, nptr = 0;
    getSymtab().getStringStats(std_str().c_str(), &ninst, &nptr);
    printf("number of shared_string instances  referencing this std::string = %li\n", ninst);
    printf("number of distinct shared_ptr sets referencing this std::string = %li\n", nptr);
}

// print debug information on the symbol table: statistics, topology of the pointer graph
void shared_string::debugSymtabStats()
{

    // topology: string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
    typedef std::map<std::string, std::map<std::shared_ptr<Symbol>*, size_t> > map_type;
    typedef std::map<std::shared_ptr<Symbol>*, size_t> imap_type;
    map_type topology;
    getSymtab().getTopology(topology);

    size_t pcnt = 0, rcnt = 0, scnt = 0;
    for (map_type::iterator it1 = topology.begin(); it1 != topology.end(); ++it1)
    {

        scnt++;
        pcnt = 0; rcnt = 0;

        printf("string entry    %li \"%s\"", scnt, it1->first.c_str());

        for (imap_type::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
            pcnt++;
            rcnt += it2->second;
            printf("    shared_ptr set nr. = %li, address = %li, reference count = %li\n", pcnt, it2->second, it2->second);
        }
        printf("total reference count = %li\n", rcnt);
    }

    printf("Total memory occupied by symbol table: %li byte\n", getSymtab().getTotalMem());

}


}      // end namespace gmx
