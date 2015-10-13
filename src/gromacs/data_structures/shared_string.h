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

#ifndef gmx_shared_string_h
#define gmx_shared_string_h

#include <string>
#include <cstring>
#include <map>
#include <iostream>
#include <ios>
#include <stdexcept>

#include "gromacs/data_structures/symtab.h"

namespace gmx
{

/*!
   \class shared_string shared_string.h "gromacs/data_structures/shared_string.h"

   \brief
    non-redundant storage of strings with std::string-like usage,
    enables global replacement

   \param[in] char*           pointer to a char string to be stored or looked up.
   \param[in] std::string     a string to be stored or looked up.
   \param[in] shared_string   a std::string to be stored or looked up.

   \returns shared_string upon inserting, retrieving, modifying an entry

   \throws std::bad_alloc if out of memory

   \throws std::logic_error in case of an (obviously) broken implementation

   A C++ replacement for the C implementation of the symbol table "t_symtab"
   for convenient use just like a std::string and automatic, clean memory
   management through the use of std::shared_ptr.

   The actual string is stored non-redundantly in the underlying symbol table
   t_symtab, shared among all instances of shared_string that represent the same
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

   \ingroup module_data_structures
   \inpublicapi
 */
class shared_string
{
    public:
// ugly, but since I need a common typedef and there are no template typedefs below C++11
#if  __cplusplus >= 201103L
 #ifdef gmx_symtab_cpp11_h
        //! basic shared_pointer type in presence of C++11 support, using the C++11 version of the symbol table
        template <class T>   using shared_ptr_type =   std::shared_ptr<T>;
 #else
        //! basic shared_pointer type in presence of C++11 support, but using the C++98 version of the symbol table
        template <class T>   using shared_ptr_type = boost::shared_ptr<T>;
 #endif
        //! shared_pointer type used within shared_string to reference the string (in presence of C++11 support)
        typedef shared_ptr_type<shared_ptr_type<gmx::t_symbol> >     str_ptr_type;
#else
        //! shared_pointer type used within shared_string to reference the string (in absence of C++11 support)
        typedef boost::shared_ptr<boost::shared_ptr<gmx::t_symbol> > str_ptr_type;
#endif
        //! the default constructor
        shared_string();
        //! copy constructor
        //! \param[in]   a   template shared_string to be copied
        shared_string(const shared_string &a);
        //! constructor from a string represented by std::string
        //! \param[in]   s   template string
        shared_string(const std::string &s);
        //! constructor from a string represented by char*
        //! \param[in]   s   template string
        shared_string(const char * const s);
        //! the destructor
        ~shared_string()
        {}
        //! replace content of this shared_string and all other instances referring to the same string with the content of *sptr
        void replace_all(const shared_string * const sptr) { replace_all(sptr->std_str()); }
        //! replace content of this shared_string and all other instances referring to the same string with the content of s
        void replace_all(const shared_string &s)
        {
            replace_all(s.std_str());
        }
        //! replace content of this shared_string and all other instances referring to the same string with the content of *s
        void replace_all(const std::string * const s) {replace_all(*s); }
        //! replace content of this shared_string and all other instances referring to the same string with the content of csptr
        //! \param[in]   csptr   pointer to the replacement string
        void replace_all(const char * const csptr)
        {
            if (csptr != nullptr)
            {
                return replace_all(std::string(csptr));
            }
            else
            {
                return;
            }
        }
        //! replace content of this shared_string and all other instances referring to the same string with the content of s
        //! \param[in]   s   the replacement string
        void replace_all(const std::string &s);
        //! defragment the symbol table, remove unused entries and unused buffers
        void clean();
        //! explicitly convert shared_string to std::string
        std::string std_string() const
        {
            if (str_ptr != nullptr && *str_ptr != nullptr)
            {
                return std::string((*str_ptr)->get_ptr());
            }
            else
            {
                return std::string();
            }
        }
        //! explictly convert shared_string to std::string
        std::string std_str() const { return std_string(); }
        //! provide this one as implicit conversion operator for convenient interchangeable use of std::string and shared_string
        operator std::string() const { return std_string(); }
        //! explicitly convert shared_string to const char* / (C-string)
        const char * c_str() const
        {
            if (str_ptr != nullptr)
            {
                return (*str_ptr)->get_ptr();
            }
            else
            {
                return nullptr;
            }
        }

        //! as for std::string, no implicit conversion to char*, usable via static_cast<const char*>(shared_string);
#if  __cplusplus >= 201103L
        explicit
#endif
                 operator const char*() const { return c_str(); }

        // begin of operators and operator related functions ...
        //! assignment operator, does not affect other shared_string instances
        shared_string &operator=(const shared_string &rhs)
        {
            if (this != &rhs)
            {
                // free old resources if not referenced by another shared_string instance
                // and reference rhs resource instead
                str_ptr = rhs.str_ptr;
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
            return (std_string() == rhs.std_string());
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
            return (std_string() < rhs.std_string());
        }
        //! comparison operator gmx::shared_string > gmx::shared_string, uses the corresponding operator of std::string
        //! \param[in]   rhs   right-hand side of the comparison
        inline bool operator>(const shared_string &rhs) const
        {
            return (std_string() > rhs.std_string());
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
            return shared_string(std::string((*str_ptr)->get_ptr()) + std::string((*rhs.str_ptr)->get_ptr()));
        }
        //! concatenate std::string and gmx::shared_string
        //! \param[in]   rhs   string to be concatenated to the end of this shared_string
        shared_string operator+(const std::string &rhs)
        {
            return shared_string(std::string((*str_ptr)->get_ptr()) + rhs);
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
#if  __cplusplus >= 201103L
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
#endif
        // stream operators
        //! ostream operator, inserts shared_string into an ostream
        friend std::ostream &operator<<(std::ostream &os, const shared_string &s);
        //! istream operator to read a shared_string from an istream
        friend std::istream &operator>>(std::istream &is, shared_string &s);
        // ... end of operators and operator related functions
        //! total number of shared_string instances sharing the same string content
        size_t get_use_count()
        {
            size_t ninst = 0;
            size_t nptr  = 0;
            get_string_stats(&ninst, &nptr);
            return ninst;
        }
        //! debugging function that returns the number of shared_ptr_type<t_symbol> sets in use by the symbol table
        size_t get_ptr_count()
        {
            size_t ninst = 0;
            size_t nptr  = 0;
            get_string_stats(&ninst, &nptr);
            return nptr;
        }
        /*! \brief debugging function that displays information on the global usage of the string content of *this

            \param[out] *ninst   number of shared_string instances referencing the same string content as *this
            \param[out] *nptr    number of distinct shared_string sets referencing the same string content as *this
                                (can be >1 after renaming) if the final string was already present in the symbol table
        */
        void get_string_stats(size_t *ninst, size_t *nptr);

        //! string length excluding the terminal '\0'
        size_t length() const { return string_length((*str_ptr)->get_ptr()) - 1; }

        //! total number of shared_string instances sharing the same outer shared_ptr<shared_ptr<t_symbol> > set
        size_t get_use_count_p() const { return str_ptr.use_count(); }

        //! total number of outer shared_ptr<shared_ptr<t_symbol> > sets sharing the same inner shared_ptr<t_symbol> in the symbol table
        size_t get_use_count_pp() const
        {
            if (str_ptr != nullptr)
            {
                return str_ptr->use_count();
            }
            else
            {
                return (size_t)0;
            }
        }

        //! get the amount of memory currently occupied by the symbol table (in bytes)
        size_t get_total_mem();

        //! get the number of strings currently stored in the symbol table
        size_t get_nr();

        //! get the number of buffers currently used by the symbol table
        size_t get_nbuf();

        /*! \brief free the symbol table, leaving the entries themselves intact

            will seldom make much sense because most of the memory will
            be occupied by the shared_string instances/share_pointers themselves */
        void free_symtab();

        //! print debug information regarding this shared_string and other shared_string instances referencing the same string content
        void debug_stats();

        //! print debug information on the symbol table: statistics, topology of the pointer graph
        void debug_symtab_stats();
    private:
        // private member functions
        //! read from istream
        inline std::istream &read(std::istream &ist);
        //! insert shared_string into ostream
        inline std::ostream &print(std::ostream &os) const
        {
            if (str_ptr != nullptr && *str_ptr != nullptr)
            {
                const char * const c_str = (*str_ptr)->get_ptr();
                return os << c_str;
            }
            else
            {
                return os;
            }
        }
        // private member variables
        /*! \brief a pointer to the central pointer in the symbol table that references the string

            shared_ptr<shared_ptr<t_symbol> > str_ptr (and str_ptr of other equivalent shared_string instances)
             -> shared_ptr<t_symbol> within symtab (and possibly other, non-equivalent sets shared_ptr<t_symbol> if renaming took place)
              -> t_symbol storing the actual string (as array of chars), and the number of non-equivalent sets shared_ptr<t_symbol> referencing the symbol
         */
        shared_ptr_type<shared_ptr_type<gmx::t_symbol> > str_ptr;
        /*! \brief
             the symbol table,
             Linking, in the header-only case, requires this wrapper function
             C++11 guarantees thread-safety with a single, lazy initialization.
        */
        inline t_symtab &get_symtab() { static t_symtab static_symtab; return static_symtab; }
};



/*! \brief comparison operators, for '>' and '<', use the corresponding operators of std::string

    operators for shared_string on the left-hand-side and std::string on right-hand-side
    are apparently implicitly converted through the assignment operator of shared_string */
/*!  std::string    comparison operator    shared_string */
//! comparison operator std::string == gmx::shared_string, uses the corresponding operator of std::string
inline bool operator==(const std::string &lhs, const shared_string &rhs)
{
    return (lhs == rhs.std_string());
}
//! comparison operator std::string != gmx::shared_string, uses the corresponding operator of std::string
inline bool operator!=(const std::string &lhs, const shared_string &rhs)
{
    return !(lhs == rhs.std_string());
}
//! comparison operator std::string < gmx::shared_string, uses the corresponding operator of std::string
inline bool operator<(const std::string &lhs, const shared_string &rhs)
{
    return (lhs < rhs.std_string());
}
//! comparison operator std::string > gmx::shared_string, uses the corresponding operator of std::string
inline bool operator>(const std::string &lhs, const shared_string &rhs)
{
    return (lhs > rhs.std_string());
}
//! comparison operator std::string >= gmx::shared_string, uses the corresponding operator of std::string
inline bool operator>=(const std::string &lhs, const shared_string &rhs)
{
    return !(lhs < rhs.std_string());
}
//! comparison operator std::string <= gmx::shared_string, uses the corresponding operator of std::string
inline bool operator<=(const std::string &lhs, const shared_string &rhs)
{
    return !(lhs > rhs.std_string());
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
    return str + sstr.std_string();
}
//! construct and assign std::string form concatenation of std::string and shared_string
//! \param[in]   str    first part  of the concatenation
//! \param[in]   sstr   second part of the concatenation
inline std::string &operator+=(std::string &str, const shared_string &sstr)
{
    return str += sstr.std_string();
}

/*! \brief global replacement of all shared_strings whose content is equal to
    before with after via the symtab, argument types T1, T2 must be implicitly
    convertible to shared_string

    \param[in] before  string to match
    \param[in] after   replacement string
*/
template<typename T1, typename T2> void replace_all_shared_strings(T1 &before, T2 &after)
{
    shared_string bstr = before;
    shared_string astr = after;
    bstr.replace_all(astr.std_string());
}

} // end namespace gmx

#endif // end headers
