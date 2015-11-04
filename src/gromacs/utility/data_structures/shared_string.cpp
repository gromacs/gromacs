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

// replace content and group ID of this shared_string and all other instances within the
// same string group and referring to the same string with the content of s and
void shared_string::replaceAll(const std::string &s, const long int group)
{
    // make sure that all entries with group IDs < 0 are added to the same default group
    long int gid_after = group;
    if (gid_after < 0)
    {
        gid_after = -1;
    }

    if (strPtr_ == nullptr)
    {
        strPtr_ = getSymtab().nonconstRenameEntry("", (long int)-1, s.c_str(), gid_after);
    }
    else
    {
        char         * tmpp               = (*strPtr_)->getPtr();
        const long int gid_before         = (*strPtr_)->getGroup();

        if (tmpp == nullptr)
        {
            strPtr_ = getSymtab().nonconstRenameEntry("", gid_before, s.c_str(), gid_after);
        }
        else
        {
            // reassigning the shared_ptr can save memory if *this is the
            // only shared_string instance referencing the current string
            strPtr_ = getSymtab().nonconstRenameEntry(tmpp, gid_before, s.c_str(), gid_after);
        }
    }
}

// append contents of gmx::shared_string rhs to shared_string
shared_string &shared_string::operator+=(const shared_string &rhs)
{
    std::string    tmpstr(std::string((*strPtr_)->getPtr()) + std::string((*rhs.strPtr_)->getPtr()));
    const long int gid = (*strPtr_)->getGroup();
    strPtr_.reset();
    strPtr_ = getSymtab().nonconstInsertEntry(tmpstr, gid);
    return *this;
}

// append contents of std::string rhs to shared_string
shared_string &shared_string::operator+=(const std::string &rhs)
{
    // by not directly calling insert_entry, it is ensured that the old entry position
    // is already free for the new version if *this was the only shared_ptr instance
    // referencing the old entry
    std::string    tmpstr(std::string((*strPtr_)->getPtr()) + rhs);
    const long int gid = (*strPtr_)->getGroup();
    strPtr_.reset();
    strPtr_ = getSymtab().nonconstInsertEntry(tmpstr, gid);
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
    long int gid = -1;
    // keep the current string group if there is one assigned already
    if (strPtr_ != nullptr)
    {
        (*strPtr_)->getGroup();
    }
    strPtr_ = getSymtab().nonconstInsertEntry(rhs, gid);
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
    getSymtab().getStringStats(std_str().c_str(), getGroup(), ninst, nptr);
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
    getSymtab().getStringStats(std_str().c_str(), getGroup(), &ninst, &nptr);
    printf("number of shared_string instances  referencing this std::string = %lu\n", ninst);
    printf("number of distinct shared_ptr sets referencing this std::string = %lu\n", nptr);
}

// print debug information on the symbol table: statistics, topology of the pointer graph
void shared_string::debugSymtabStats()
{

    // topology: string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
    typedef std::shared_ptr<Symbol>*         optr_type;
    typedef std::map<optr_type, size_t>      imap_type;
    typedef std::map<std::string, imap_type> mmap_type;
    typedef std::map<long int, mmap_type>    omap_type;
    omap_type topology;
    getSymtab().getTopology(topology);

    for (omap_type::iterator it1 = topology.begin(); it1 != topology.end(); ++it1)
    {
        printf("\n--------------------\n");
        printf("string group %li:\n", it1->first);

        size_t pcnt = 0, rcnt = 0, scnt = 0;
        for (mmap_type::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {

            scnt++;
            pcnt = 0; rcnt = 0;

            printf("string entry    %lu \"%s\"", scnt, it2->first.c_str());

            for (imap_type::iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
            {
                pcnt++;
                rcnt += it3->second;
                printf("    shared_ptr set nr. = %li, address = %p, reference count = %li\n", pcnt, static_cast<void*>(it3->first), it3->second);
            }
            printf("total reference count = %lu\n", rcnt);
        }
    }

    printf("Total memory occupied by symbol table: %lu byte\n", getSymtab().getTotalMem());

}


}      // end namespace gmx
