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
    \inpublicapi

   \brief
    Class shared_string stores std::strings non-redundantly (thread-safe) and trimmed
    to the chars actually used in a symbol table. shared_string can be used much like
    std::string plus allowing for a global replacement of a std::string, e.g., for
    renaming all residues of a certain type in GROMACS.

    Identical strings of characters are stored non-redundantly in the symbol table and
    shared among all instances of shared_string that represent the same string via
    shared_pointers. This header selects between the C++11 and C++98 variants of the
    symbol table according to the availability of the C++11 standard as signalled by
    the compiler via preprocessor macros. The C++11 version is closer to the aim of
    automatic memory management in letting string entries that are not used outside
    the symbol table anymore also automatically vanish from the symbol table.

   \author
    R. Thomas Ullmann <tullman@gwdg.de>

   \copyright
    GROMACS license

   \date Mar 2015
 */

#ifndef GMX_UTILITY_SHARED_STRING_H
#define GMX_UTILITY_SHARED_STRING_H

#include <cstring>

#include <ios>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "gromacs/utility/data_structures/symtab_cpp.h"


namespace gmx
{

/*!
   \class shared_string

   \brief
    non-redundant storage of strings with std::string-like usage,
    enables global replacement

   \param[in] char*           pointer to a char string to be stored or looked up.
   \param[in] std::string     a string to be stored or looked up.
   \param[in] shared_string   a std::string to be stored or looked up.

   \returns shared_string upon inserting, retrieving, modifying an entry

   \throws std::bad_alloc if out of memory

   \throws std::logic_error in case of an (obviously) broken implementation

   A C++ replacement for the C implementation of the symbol table "Symtab"
   for convenient use just like a std::string and automatic, clean memory
   management through the use of std::shared_ptr.

   The actual string is stored non-redundantly in the underlying symbol table
   Symtab, shared among all instances of shared_string that represent the same
   string of characters.

   Shared_string provides some of the functionality of std::string, while
   being memory saving and intrinsically thread save -- possibly at the
   expense of some performance in write operations to the symbol table.

   Accesses to the symbol table are currently not lock-free for simplicity and
   because handling the std::strings will most likely not be performance critical
   in my applications. Simple usage of the std::string content without altering it
   should not be affected because the symbol table is not involved or affected.

   shared_string is intended for storing a limited number of unique entries/strings
   like the names of residues, atoms, sites and site forms in GROMACS.

   \ingroup module_utility
   \inpublicapi
 */
class shared_string
{
    public:
        //! shared_pointer type used within shared_string to reference the string
        typedef std::shared_ptr<std::shared_ptr<gmx::Symbol> > strPtr_type;

        //! the default constructor
        shared_string()
        {
            strPtr_ = getSymtab().nonconstInsertEntry("");
        }

        //! copy constructor
        //! \param[in]   a   template shared_string to be copied
        shared_string(const shared_string &a) : strPtr_(a.strPtr)
        {
        }

        //! constructor from a string represented by std::string
        //! \param[in]   s   template string
        shared_string(const std::string &s)
        {
            strPtr_ = getSymtab().nonconstInsertEntry(s.c_str());
        }

        //! constructor from a string represented by std::string
        //! \param[in]   s       template string
        //! \param[in]   group   group to add the shared_string to
        shared_string(const std::string &s, const long int group)
        {
            strPtr_ = getSymtab().nonconstInsertEntry(s.c_str(), group);
        }

        //! constructor from a string represented by char*
        //! \param[in]   s   template string
        shared_string(const char * const s)
        {
            if (s != nullptr)
            {
                strPtr_ = getSymtab().nonconstInsertEntry(s);
            }
            else
            {
                strPtr_ = getSymtab().nonconstInsertEntry("");
            }
        }

        //! constructor from a string represented by char*
        //! \param[in]   s       template string
        //! \param[in]   group   group to add the shared_string to
        shared_string(const char * const s, const long int group)
        {
            if (s != nullptr)
            {
                strPtr_ = getSymtab().nonconstInsertEntry(s, group);
            }
            else
            {
                strPtr_ = getSymtab().nonconstInsertEntry("", group);
            }
        }

        //! the destructor
        ~shared_string() {}

        //! replace content of this shared_string and all other instances referring to the same string with the content of *sptr
        //! set the default group ID equal to -1 if none is assigned yet, otherwise keep the current group ID
        //! \param[in]   sptr   pointer to the replacement string
        void replaceAll(const shared_string * const sptr)
        { replaceAll(sptr->std_str(), getGroup()); }
        //! replace content of this shared_string and all other instances referring to the same string with the content of s
        //! set the default group ID equal to -1 if none is assigned yet, otherwise keep the current group ID
        //! \param[in]   s   the replacement string
        void replaceAll(const shared_string &s)
        {
            replaceAll(s.std_str(), getGroup());
        }
        //! replace content of this shared_string and all other instances referring to the same string with the content of *s
        //! set the default group ID equal to -1 if none is assigned yet, otherwise keep the current group ID
        //! \param[in]   s   the replacement string
        void replaceAll(const std::string * const s)
        { replaceAll(*s); }
        //! replace content of this shared_string and all other instances referring to the same string with the content of csptr
        //! set the default group ID equal to -1 if none is assigned yet, otherwise keep the current group ID
        //! \param[in]   csptr   pointer to the replacement string
        void replaceAll(const char * const csptr)
        {
            if (csptr != nullptr)
            {
                return replaceAll(std::string(csptr), getGroup());
            }
            else
            {
                return;
            }
        }
        //! replace content of this shared_string and all other instances referring to the same string with the content of s
        //! set the default group ID equal to -1 if none is assigned yet, otherwise keep the current group ID
        //! \param[in]   s   the replacement string
        void replaceAll(const std::string &s)
        {
            return replaceAll(s, getGroup());
        }

        //! replace content of this shared_string and all other instances referring to the same string with the content of *sptr
        //! \param[in]   sptr   pointer to the replacement string
        //! \param[in]   group_after   the replacement group ID
        void replaceAll(const shared_string * const sptr, long int group_after)
        { replaceAll(sptr->std_str(), group_after); }
        //! replace content of this shared_string and all other instances referring to the same string with the content of s
        //! \param[in]   s             the replacement string
        //! \param[in]   group_after   the replacement group ID
        void replaceAll(const shared_string &s, long int group_after)
        {
            replaceAll(s.std_str(), group_after);
        }
        //! replace content of this shared_string and all other instances referring to the same string with the content of *s
        //! \param[in]   s             the replacement string
        //! \param[in]   group_after   the replacement group ID
        void replaceAll(const std::string * const s, long int group_after)
        { replaceAll(*s, group_after); }
        //! replace content of this shared_string and all other instances referring to the same string with the content of csptr
        //! \param[in]   csptr         pointer to the replacement string
        //! \param[in]   group_after   the replacement group ID
        void replaceAll(const char * const csptr, long int group_after)
        {
            if (csptr != nullptr)
            {
                return replaceAll(std::string(csptr), group_after);
            }
            else
            {
                return;
            }
        }
        //! replace content of this shared_string and all other instances referring to the same string with the content of s
        //! \param[in]   s             the replacement string
        //! \param[in]   group_after   the replacement group ID
        void replaceAll(const std::string &s, long int group_after);

        //! defragment the symbol table, remove unused entries and unused buffers
        void clean();

        //! explicitly convert shared_string to std::string
        //! alias of std_str() in GROMACS naming convention
        std::string stdString() const
        {
            if (strPtr_ != nullptr && *strPtr_ != nullptr)
            {
                return std::string((*strPtr_)->getPtr());
            }
            else
            {
                return std::string();
            }
        }
        //! explictly convert shared_string to std::string
        //! this function is named std_str as in std::string instead of stdStr for consistency with c_str()
        std::string std_str() const { return stdString(); }
        //! provide this one as implicit conversion operator for convenient interchangeable use of std::string and shared_string
        operator std::string() const { return stdString(); }

        //! explicitly convert shared_string to const char* / (C-string)
        //! this function is named c_str as in std::string instead of cStr to avoid confusion
        const char * c_str() const
        {
            if (strPtr_ != nullptr)
            {
                return (*strPtr_)->getPtr();
            }
            else
            {
                return nullptr;
            }
        }
        //! explicitly convert shared_string to const char* / (C-string)
        //! alias of c_str() in GROMACS naming convention
        const char * cString() const { return c_str(); }

        //! as for std::string, but no implicit conversion to char*, usable via static_cast<const char*>(shared_string);
        explicit operator const char*() const { return c_str(); }

        //! assignment operator, does not affect other shared_string instances
        shared_string &operator=(const shared_string &rhs)
        {
            if (this != &rhs)
            {
                // free old resources if not referenced by another shared_string instance
                // and reference rhs resource instead
                strPtr_ = rhs.strPtr_;
            }
            return *this;
        }

        //! assign left-hand-side shared_string individually to right-hand-side string
        //! \param[in]   rhs   string content to be assigned to this shared_string
        shared_string &operator=(const std::string &rhs);

        //! assign left-hand-side shared_string individually to right-hand-side shared string
        //! \param[in]   rhs   string content to be assigned to this shared_string
        shared_string &operator=(const char * const rhs);

        //! read a string from an istream into the shared_string -- shared_string << istr
        //! \param[in]   ist   istream from which string content is to be read
        std::istream &operator<<(std::istream &ist);

        // comparison operators, for '>' and '<', use the corresponding operators of std::string
        // shared_string on right-hand-side
        //! comparison operator gmx::shared_string == gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator==(const shared_string &rhs) const
        {
            return (stdString() == rhs.stdString());
        }
        //! comparison operator gmx::shared_string != gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator!=(const shared_string &rhs) const
        {
            return !(*this == rhs);
        }
        //! comparison operator gmx::shared_string < gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator<(const shared_string &rhs) const
        {
            return (stdString() < rhs.stdString());
        }
        //! comparison operator gmx::shared_string > gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator>(const shared_string &rhs) const
        {
            return (stdString() > rhs.stdString());
        }
        //! comparison operator gmx::shared_string >= gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator>=(const shared_string &rhs) const
        {
            return !(*this < rhs);
        }
        //! comparison operator gmx::shared_string <= gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator<=(const shared_string &rhs) const
        {
            return !(*this > rhs);
        }

        // string manipulation operators
        //! concatenate
        //! \param[in]   rhs   string to be concatenated to the end of this shared_string
        shared_string operator+(const shared_string &rhs)
        {
            return shared_string(std::string((*strPtr_)->getPtr()) + std::string((*rhs.strPtr_)->getPtr()));
        }
        //! concatenate std::string and gmx::shared_string
        //! \param[in]   rhs   string to be concatenated to the end of this shared_string
        shared_string operator+(const std::string &rhs)
        {
            return shared_string(std::string((*strPtr_)->getPtr()) + rhs);
        }

        // modification and assignment
        // by not directly calling insert_entry, it is ensured that the old entry position
        // is already free for the new version if *this was the only shared_ptr instance
        // referencing the old entry
        //! append contents of gmx::shared_string rhs to shared_string
        //! \param[in]   rhs   string to be appended to the content of this shared_string
        shared_string &operator+=(const shared_string &rhs);
        //! append contents of std::string rhs to shared_string
        //! \param[in]   rhs   string to be appended to the content of this shared_string
        shared_string &operator+=(const std::string &rhs);
        //! in principle possible, but potentially complicated/ambiguous, needed?
        shared_string operator-() = delete;
        // explicitly declare increment and decrement operators deleted
        // post-increment
        //! make clear that shared_str and child classes don't offer post-increment &operator--
        shared_string &operator++() = delete;
        //! make clear that shared_str and child classes don't offer post-decrement &operator--
        shared_string &operator--() = delete;
        // pre-increment
        //! make clear that shared_str and child classes don't offer pre-increment operator++
        shared_string operator++(int) = delete;
        //! make clear that shared_str and child classes don't offer pre-decrement operator--
        shared_string operator--(int) = delete;

        // stream operators
        //! ostream operator, inserts shared_string into an ostream
        friend std::ostream &operator<<(std::ostream &os, const shared_string &s);
        //! istream operator to read a shared_string from an istream
        friend std::istream &operator>>(std::istream &is, shared_string &s);

        // ... end of operators and operator related functions

        //! retrieve the ID of the string group of which this shared_string is member
        long int getGroup()
        {
            // return the default group ID if none is assigned yet
            long int gid = -1;
            if (strPtr_ != nullptr)
            {
                gid = (*strPtr_)->getGroup();
            }
            return gid;
        }

        //! set the string group ID to change the string group membership
        //! \param[in]   group   new group ID
        void setGroup(const long int group)
        {
            // set default group ID -1 for any input group ID < 0
            long int new_gid = group;
            if (new_gid < 0)
            {
                new_gid = -1;
            }
            // Can not simply set the group ID of the referenced Symbol because this
            // would affect other instances referencing the same Symbol, too.
            strPtr_ = getSymtab().nonconstInsertEntry(std_str(), new_gid);
        }

        //! set the string group ID to change the string group membership, return the previous group ID
        //! \param[in]   group   new group ID
        long int exchangeGroup(const long int group)
        {
            long int old_gid = getGroup();
            // set default group ID -1 for any input group ID < 0
            long int new_gid = group;
            if (new_gid < 0)
            {
                new_gid = -1;
            }
            // Can not simply set the group ID of the referenced Symbol because this
            // would affect other instances referencing the same Symbol, too.
            strPtr_ = getSymtab().nonconstInsertEntry(std_str(), new_gid);
            return old_gid;
        }

        //! total number of shared_string instances sharing the same string content
        size_t getUseCount()
        {
            size_t ninst = 0;
            size_t nptr  = 0;
            getStringStats(&ninst, &nptr);
            return ninst;
        }
        //! debugging function that returns the number of shared_ptr_type<Symbol> sets in use by the symbol table
        size_t getPtrCount()
        {
            size_t ninst = 0;
            size_t nptr  = 0;
            getStringStats(&ninst, &nptr);
            return nptr;
        }
        //! \brief debugging function that displays information on the global usage of the string content of *this
        //!
        //! \param[out] *ninst   number of shared_string instances referencing the same string content as *this
        //! \param[out] *nptr    number of distinct shared_string sets referencing the same string content as *this
        //!                      (can be >1 after renaming) if the final string was already present in the symbol table
        void getStringStats(size_t *ninst, size_t *nptr);

        //! string length excluding the terminal '\0'
        size_t length() const { return std::strlen((*strPtr_)->getPtr()); }

        //! total number of shared_string instances sharing the same outer shared_ptr<shared_ptr<Symbol> > set
        size_t getUseCountP() const { return strPtr_.use_count(); }
        //! alias of getUseCountP for consistency with shared_ptr::use_count()
        size_t use_count() const { return getUseCountP(); }

        //! total number of outer shared_ptr<shared_ptr<Symbol> > sets sharing the same inner shared_ptr<Symbol> in the symbol table
        size_t getUseCountPP() const
        {
            if (strPtr_ != nullptr)
            {
                return strPtr_->use_count();
            }
            else
            {
                return (size_t)0;
            }
        }

        //! get the amount of memory currently occupied by the symbol table (in bytes)
        size_t getTotalMem();

        //! get the number of strings currently stored in the symbol table
        size_t getNr();

        //! get the number of buffers currently used by the symbol table
        size_t getNbuf();

        /*! \brief free the symbol table, leaving the entries themselves intact

            will seldom make much sense because most of the memory will
            be occupied by the shared_string instances/share_pointers themselves */
        void freeSymtab();

        //! print debug information regarding this shared_string and other shared_string instances referencing the same string content
        void debugStats();

        //! print debug information on the symbol table: statistics, topology of the pointer graph
        void debugSymtabStats();
    private:
        // private member functions
        //! read from istream
        inline std::istream &read(std::istream &ist);
        //! insert shared_string into ostream
        inline std::ostream &print(std::ostream &os) const
        {
            if (strPtr_ != nullptr && *strPtr_ != nullptr)
            {
                const char * const c_str = (*strPtr_)->getPtr();
                return os << c_str;
            }
            else
            {
                return os;
            }
        }
        // private member variables
        /*! \brief a pointer to the central pointer in the symbol table that references the string

            shared_ptr<shared_ptr<Symbol> > strPtr_ (and strPtr_ of other equivalent shared_string instances)
             -> shared_ptr<Symbol> within symtab (and possibly other, non-equivalent sets shared_ptr<Symbol> if renaming took place)
              -> Symbol storing the actual string (as array of chars), and the number of non-equivalent sets shared_ptr<Symbol> referencing the symbol
         */
        strPtr_type strPtr_;
        /*! \brief
             the symbol table,
             Linking, in the header-only case, requires this wrapper function
             C++11 guarantees thread-safety with a single, lazy initialization.
         */
        inline Symtab &getSymtab() { static Symtab static_symtab; return static_symtab; }
};

// read from istream
inline std::istream &shared_string::read(std::istream &ist)
{
    std::string tmp_str;
    ist >> tmp_str;
    if (ist.good())
    {
        getSymtab().nonconstInsertEntry(tmp_str.c_str());
    }
    return ist;
}

// read shared_string from istream
inline std::istream &operator>>(std::istream &is, shared_string &s)
{
    return s.read(is);
}

/*! \brief comparison operators, for '>' and '<', use the corresponding operators of std::string

    operators for shared_string on the left-hand-side and std::string on the right-hand-side
    are apparently implicitly converted through the assignment operator of shared_string */
/*!  std::string    comparison operator    shared_string */
//! comparison operator std::string == gmx::shared_string, uses the corresponding operator of std::string
inline bool operator==(const std::string &lhs, const shared_string &rhs)
{
    return (lhs == rhs.stdString());
}
//! comparison operator std::string != gmx::shared_string, uses the corresponding operator of std::string
inline bool operator!=(const std::string &lhs, const shared_string &rhs)
{
    return !(lhs == rhs.stdString());
}
//! comparison operator std::string < gmx::shared_string, uses the corresponding operator of std::string
inline bool operator<(const std::string &lhs, const shared_string &rhs)
{
    return (lhs < rhs.stdString());
}
//! comparison operator std::string > gmx::shared_string, uses the corresponding operator of std::string
inline bool operator>(const std::string &lhs, const shared_string &rhs)
{
    return (lhs > rhs.stdString());
}
//! comparison operator std::string >= gmx::shared_string, uses the corresponding operator of std::string
inline bool operator>=(const std::string &lhs, const shared_string &rhs)
{
    return !(lhs < rhs.stdString());
}
//! comparison operator std::string <= gmx::shared_string, uses the corresponding operator of std::string
inline bool operator<=(const std::string &lhs, const shared_string &rhs)
{
    return !(lhs > rhs.stdString());
}

//! write shared_string to ostream -- ostream << shared_string
//! \param[in]   os   target ostream for output
//! \param[in]   s    shared_string to be inserted into os
inline std::ostream &operator<<(std::ostream &os, const shared_string &s)
{
    return s.print(os);
}

//! read shared_string from istream -- istream >> shared_string
//! \param[in]   is   source istream
//! \param[in]   s    target shared_string for the content of is
inline std::istream &operator>>(std::istream &is, shared_string &s);

//! construct std::string form concatenation of std::string and shared_string
//! \param[in]   str    first part  of the concatenation
//! \param[in]   sstr   second part of the concatenation
inline std::string operator+(const std::string &str, const shared_string &sstr)
{
    return str + sstr.stdString();
}
//! construct and assign std::string form concatenation of std::string and shared_string
//! \param[in]   str    first part  of the concatenation
//! \param[in]   sstr   second part of the concatenation
inline std::string &operator+=(std::string &str, const shared_string &sstr)
{
    return str += sstr.stdString();
}

/*! \brief global replacement of all shared_strings whose content is equal to
    before with after via the symtab, argument types T1, T2 must be implicitly
    convertible to shared_string

    \tparam    T1    data type of string before
    \tparam    T2    data type of string after

    \param[in] before  string to match
    \param[in] after   replacement string
 */
template<typename T1, typename T2> void replaceAllSharedStrings(T1 &before, T2 &after)
{
    shared_string bstr = before;
    shared_string astr = after;
    bstr.replaceAll(astr.stdString());
}

/*! \brief global replacement of all shared_strings whose content is equal to
    before with after via the symtab, argument types T1, T2 must be implicitly
    convertible to shared_string

    \tparam    T1    data type of string before
    \tparam    T2    data type of string after

    \param[in] before        string to match
    \param[in] group_before  group of before to match
    \param[in] after         replacement string
    \param[in] group_after   replacement group
 */
template<typename T1, typename T2> void replaceAllSharedStrings(T1 &before, long int group_before, T2 &after, long int group_after)
{
    // make sure that all negative group IDs are mapped to the default group -1
    long int gid_before = group_before;
    long int gid_after  = group_after;
    if (gid_before <  0)
    {
        gid_before = -1;
    }
    if (gid_after  <  0)
    {
        gid_before = -1;
    }
    shared_string bstr = before;
    bstr.setGroup(gid_before);
    shared_string astr = after;
    bstr.replaceAll(astr.stdString(), gid_after);
}

}      // end namespace gmx

#endif // end header
