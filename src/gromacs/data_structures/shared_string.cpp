/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015 by the GROMACS development team, led by
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
    \ingroup module_data_structures
    \inpublicapi

   \brief
    Implementation file needed for linking all code parts using
    shared_string to a common shared_string object.

   \author
    R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \copyright
    GROMACS license

   \date Sep 2015
 */

#include "shared_string.h"


namespace gmx
{

// the default constructor
shared_string::shared_string()
{
    str_ptr = get_symtab().nonconst_insert_entry("");
}

// copy constructor
shared_string::shared_string(const shared_string &a)
{
    str_ptr = a.str_ptr;
}

// constructor from a string represented by std::string
shared_string::shared_string(const std::string &s)
{
    str_ptr = get_symtab().nonconst_insert_entry(s.c_str());
}

// constructor from a string represented by char*
shared_string::shared_string(const char * const s)
{
    if (s != nullptr)
    {
        str_ptr = get_symtab().nonconst_insert_entry(s);
    }
    else
    {
        str_ptr = get_symtab().nonconst_insert_entry("");
    }
}

// replace content of this shared_string and all other instances referring to the same string with the content of s
void shared_string::replace_all(const std::string &s)
{
    if (str_ptr == nullptr)
    {
        str_ptr = get_symtab().nonconst_rename_entry("", s.c_str());
    }
    else
    {
        char * tmpp = (*str_ptr)->get_ptr();
        if (tmpp == nullptr)
        {
            str_ptr = get_symtab().nonconst_rename_entry("", s.c_str());
        }
        else
        {
            // reassigning the shared_ptr can save memory if *this is the
            // only shared_string instance referencing the current string
            str_ptr = get_symtab().nonconst_rename_entry(tmpp, s.c_str());
        }
    }
}

// append contents of gmx::shared_string rhs to shared_string
shared_string& shared_string::operator+=(const shared_string &rhs)
{
    std::string tmpstr(std::string((*str_ptr)->get_ptr()) + std::string((*rhs.str_ptr)->get_ptr()));
    str_ptr.reset();
    str_ptr = get_symtab().nonconst_insert_entry(tmpstr);
    return *this;
}

// append contents of std::string rhs to shared_string
shared_string& shared_string::operator+=(const std::string &rhs)
{
    // by not directly calling insert_entry, it is ensured that the old entry position
    // is already free for the new version if *this was the only shared_ptr instance
    // referencing the old entry
    std::string tmpstr(std::string((*str_ptr)->get_ptr()) + rhs);
    str_ptr.reset();
    str_ptr = get_symtab().nonconst_insert_entry(tmpstr);
    return *this;
}

// defragment the symbol table, remove unused entries and unused buffers
void shared_string::clean()
{
    return get_symtab().clean();
}

// assign left-hand-side shared_string individually to right-hand-side string
shared_string& shared_string::operator=(const std::string &rhs)
{
    str_ptr = get_symtab().nonconst_insert_entry(rhs);
    return *this;
}
// assign left-hand-side shared_string individually to right-hand-side shared string
shared_string& shared_string::operator=(const char * const rhs)
{
    str_ptr = get_symtab().nonconst_insert_entry(rhs);
    return *this;
}

// read a string from an istream into the shared_string -- shared_string << istr
std::istream& shared_string::operator<<(std::istream &ist)
{
    return read(ist);
}

// free the symbol table, leaving the entries themselves intact 
void shared_string::free_symtab()
{
    return get_symtab().free_symtab(false);
}

/* \brief debugging function that displays information on the global usage of the string content of *this

    \param[out] *ninst   number of shared_string instances referencing the same string content as *this
    \param[out] *nptr    number of distinct shared_string sets referencing the same string content as *this
                        (can be >1 after renaming) if the final string was already present in the symbol table
*/
void shared_string::get_string_stats(size_t *ninst, size_t *nptr)
{
    get_symtab().get_string_stats(std_str().c_str(), ninst, nptr);
}

// get the amount of memory currently occupied by the symbol table (in bytes)
size_t shared_string::get_total_mem()
{
    return get_symtab().get_total_mem();
}

// get the number of strings currently stored in the symbol table
size_t shared_string::get_nr()
{
    return get_symtab().get_nr();
}

// get the number of buffers currently used by the symbol table
size_t shared_string::get_nbuf()
{
    return get_symtab().get_nbuf();
}

// print debug information regarding this shared_string and other shared_string instances referencing the same string content
void shared_string::debug_stats()
{
    if (str_ptr)
    {
        printf("  1st level is OK.\n");
        if (*str_ptr)
        {
            printf("    2nd level is OK. use count %li", str_ptr.use_count());
            printf("    (*str_ptr)->get_ref()       : %li\n", (*str_ptr)->get_ref());
            printf("    (*str_ptr)->get_ptr()       : %p\n", static_cast<void*>((*str_ptr)->get_ptr()));
            if ((*str_ptr)->get_ptr())
            {
                printf("    (*str_ptr)->get_ptr() std::string: \"%s\"\n", (*str_ptr)->get_ptr());
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
    get_symtab().get_string_stats(std_str().c_str(), &ninst, &nptr);
    printf("number of shared_string instances  referencing this std::string = %li\n", ninst);
    printf("number of distinct shared_ptr sets referencing this std::string = %li\n", nptr);
}

// print debug information on the symbol table: statistics, topology of the pointer graph
void shared_string::debug_symtab_stats()
{

    // topology: string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
#ifdef gmx_symtab_cpp11_h
    typedef std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> > map_type;
    typedef std::map<std::shared_ptr<t_symbol>*, size_t> imap_type;
#else
    typedef std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> > map_type;
    typedef std::map<boost::shared_ptr<t_symbol>*, size_t> imap_type;
#endif
    map_type topology;
    get_symtab().get_topology(topology);

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

    printf("Total memory occupied by symbol table: %li byte\n", get_symtab().get_total_mem());

}

// read from istream
inline std::istream& shared_string::read(std::istream &ist)
{
    std::string tmp_str;
    ist >> tmp_str;
    if (ist.good())
    {
        get_symtab().nonconst_insert_entry(tmp_str.c_str());
    }
    return ist;
}

inline std::istream &operator>>(std::istream &is, shared_string &s)
{
    return s.read(is);
}


} // end namespace gmx
