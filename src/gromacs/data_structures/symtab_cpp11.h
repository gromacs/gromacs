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
    \ingroup module_data_structures
    \inpublicapi

   \brief
   Class t_symtab stores std::strings non-redundantly (thread-safe) and
   trimmed to the chars actually used. Global replacement is possible, e.g.,
   for renaming all residues of a certain type in GROMACS.

   The class is intended as C++ replacement for the C implementation of the
   symbol table "t_symtab" for convenient use just like a std::string and
   automatic, clean memory management.

   Only suitable for a limited number of unique entries/strings because the
   entries are stored unsorted, not allowing for time-saving lexical search.
   -> entry retrieval time grows linearly with the number of unique entries.

   The entries are not modifiable directly but only via dedicated member
   functions/operators to ensure internal consistency.

   The symbol table contains a linked list of symbol buffers t_symbuf.
   Each t_symbuf encloses a set of string entries defined by bufsize.
   Each string entry is saved as array of chars. Excessively long strings
   are truncated to a maximum length defined via max_entry_size.

   C++11 features are used for thread safety and improved memory management.
   C++11 smart pointers are used to keep track of the number of times
   entries are referenced and for automatic disposal of unreferenced entries.
   An entry ceases to exist as soon as no more external shared_ptr references
   are held.

   \author
    R. Thomas Ullmann <tullman@gwdg.de>

   \copyright
    GROMACS license

   \date Jan 2014

 */

#ifndef gmx_symtab_cpp11_h
#define gmx_symtab_cpp11_h

#include <string>
#include <cstring>
#include <map>
#include <stdexcept>

// check for C++11 compiler, should have been checked already in shared_string.h, left it here as a reminder
#if __cplusplus >= 201103L
 #include <mutex>
 #include <memory>
 #include <thread>
 #include <atomic>
#else
 #error "this symbol table variant requires a C++ compiler supporting the C++11 standard"
#endif

namespace gmx
{

//! \addtogroup module_data_structures
//! \{

/// \cond DEV

//! count the number of characters in an array of chars (including the terminating '\0')
static size_t string_length(const char * s)
{
    size_t result = 0;
    while (s != nullptr)
    {
        result++;
        if (*s == '\0')
        {
            break;
        }
        s++;
    }
    return result;
}

class t_symbol
{
    public:
        t_symbol()
        {
            ptr.store(nullptr);
            cnt.store(1);
        }
        ~t_symbol()
        {
            if (ptr != nullptr)
            {
                const char* tmp_p = ptr.exchange(nullptr);
                if (tmp_p != nullptr)
                {
                    delete [] tmp_p;
                }
            }
        }
    private:
        //! only classes constituting shared_string can access the class members
        friend class t_symbuf;
        friend class t_symtab;
        friend class shared_string;
        // private member functions:
        size_t use_count() const { return get_ref(); }
        size_t get_ref() const { return cnt.load(std::memory_order_seq_cst); }
        void   set_ref(const size_t &val) { cnt.store(val, std::memory_order_seq_cst); }
        size_t incr_ref()
        { return (cnt.fetch_add((size_t) 1, std::memory_order_seq_cst) + (size_t)1); }
        size_t decr_ref()
        { return (cnt.fetch_add((size_t)-1, std::memory_order_seq_cst) - (size_t)1); }
        char* get_ptr() { return ptr.load(std::memory_order_seq_cst); }
        void set_ptr(char* const new_ptr)
        { return ptr.store(new_ptr, std::memory_order_seq_cst); }
        char* exchange_ptr(char* const new_ptr)
        { return ptr.exchange(new_ptr, std::memory_order_seq_cst); }
        // private member variables:
        // pointer to the actual string
        std::atomic<char*>  ptr;
        // 2nd level reference count:
        // number of distinct sets of shared_ptr "spp" referencing
        // the shared pointer that holds the atomic string
        // spp_{1,...,n} -> sp -> t_symbol
        std::atomic<size_t> cnt;
};

class t_symbuf
{
    public:
        //! the default constructor
        t_symbuf() : prev()
        {
            bufsize   = 0;
            buf       = nullptr;
            next      = nullptr;
        }
        //! explicitly state that no copying is allowed
        t_symbuf(const t_symbuf &) = delete;
        //! the destructor
        ~t_symbuf()
        {
            // by default, we leave the string entries intact for further use
            // the shared_ptr will take care of their deallocation later on
            free_symbuf(false);
            next = nullptr;
        }
    private:
        //!  only t_symtab is granted direct access to the class members
        friend class t_symtab;
        // member variables
        //! maximum number of stored chars
        size_t bufsize;
        /*! pointer to the stored data (array of smart pointers to strings)
            use weak_ptr to ensure that a string is automatically destroyed
            as soon as there is no more external shared_ptr referencing it  */
        std::weak_ptr<std::shared_ptr<t_symbol> > *buf;
        //! pointer to the next t_symbuf
        std::shared_ptr<t_symbuf>                  next;
        //! pointer to the previous t_symbuf, weak_ptr to prevent cyclic blockage of deallocation
        std::weak_ptr<t_symbuf>                    prev;
        // member functions
        std::string trim_string(const char * const s, const size_t maxlen) const;
        bool match_strings(const std::string &s1, const std::string &s2) const {return match_strings(s1.c_str(), s2.c_str()); }
        bool match_strings(const char * const s1, const char * const s2) const;
        void free_symbuf(){return free_symbuf(false); }
        void free_symbuf(const bool &dealloc_entries)
        {
            if (buf != nullptr)
            {
                // on request, enforce deletion of the stored strings
                // even if there are still shared_ptrs referencing them
                if (dealloc_entries)
                {
                    for (size_t i = 0; i < bufsize; ++i)
                    {
                        const char* tmp_cptr = (buf[i].lock()->get())->exchange_ptr(nullptr);
                        delete [] tmp_cptr;
                    }
                }
                for (size_t i = 0; i < bufsize; ++i)
                {
                    // FIXME?: most likely NO!, seems to be correct that way, similar with a custom allocator
                    // without this reset, the weak_ptrs are not deallocated properly, heaven knows why
                    // not a problem with a bare weak_ptr or a single weak_ptr allocated on buf via new
                    // still necessary even if expired() == true here, why the hell is that so ...!?
                    // valgrind reports a line where allocator and deleter are saved inside the shared_ptr
                    // but that is long gone by now and should anyhow not be affected by whatever the
                    // weak_ptr does
                    buf[i].reset();
                }
                ::operator delete (buf);
                buf = nullptr;
            }
            bufsize = 0;
        }
        void alloc_symbuf(size_t * const mem)
        {
            try // allocate memory for the weak_ptrs
            {
                const size_t this_array_size = sizeof(std::weak_ptr<std::shared_ptr<t_symbol> >) * bufsize;
                buf = static_cast< std::weak_ptr<std::shared_ptr<t_symbol> >* > (::operator new(this_array_size));
                for (size_t i = 0; i < bufsize; ++i)
                {
                    // create the weak_ptrs in the allocated memory region
                    ::new (&buf[i]) std::weak_ptr<std::shared_ptr<std::shared_ptr<t_symbol> > >();
                    if (!buf[i].expired())
                    {
                        buf[i].reset();
                    }
                }
                *mem = sizeof(t_symbuf) + this_array_size;
            }
            catch (std::bad_alloc &)
            {
                fprintf (stderr, "Error in %s::alloc_symbuf():",  typeid(*this).name());
                fprintf (stderr, "Insufficient memory for allocating new symbol table buffer!\n");
                throw;
            }
            ;
        }
        // check existence of a string entry
        bool contains_entry(const char * const s) const {return (get_entry(s) != nullptr); }
        //! retrieve a string entry
        std::shared_ptr<std::shared_ptr<t_symbol> > get_entry(const char * const s) const;
        //! store a string entry
        bool insert_entry(std::shared_ptr<std::shared_ptr<t_symbol> > in);
        //! get number of strings stored and an estimate of the total amount of memory occupied by this t_symbuf
        void get_stats(unsigned int * const nr, size_t * const mem);
        void get_topology(std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> > &topology) const;
};

class t_symtab
{
    public:
        //! the default constructor
        t_symtab() : bufsize(81), max_entry_size(1024), nbuf(0), nr(0), symbuf(), mem(sizeof(*this)),
                     n_orphan_entries(0)
        {
        }
        //! a constructor that allows custom values for buffer size and maximum entry length
        t_symtab(const size_t &bs, const size_t &mes) : bufsize(bs), max_entry_size(mes), nbuf(0), nr(0), symbuf(), mem(sizeof(*this)),
                                                        n_orphan_entries(0)
        {
        }
        //! explicitly state that no copying is allowed
        t_symtab(const t_symtab &) = delete;
        //! the destructor
        ~t_symtab()
        {
            // should all happen by itself thanks to the shared_pointer linkage of the t_symbufs
            // free_symtab();
            // symbuf.reset();
        }
        //! by default, leave the strings intact, shared_ptrs will free them on destruction
        void free_symtab(){return free_symtab(false); }
        // (dealloc_entries == false) -> strings managed by the handed
        // out shared_ptrs stay intact
        void free_symtab (const bool &dealloc_entries)
        {
            std::shared_ptr<t_symbuf> tmp_buf(symbuf);
            // go along the linked list of symbufs deleting all entries,
            // or only the weak references to them if !dealloc_entries
            while (tmp_buf != nullptr)
            {
                tmp_buf->free_symbuf(dealloc_entries);
                tmp_buf = tmp_buf->next;
            }
            // initiate the successive calling of dtors along the linked list of t_symbufs
            symbuf.reset();
        }
        /*! defragment t_symbufs by moving entries towards the beginning of the linked list.
            remove left over empty t_symbuf structures */
        void clean();
        //! maximum number of elements stored per symbol buffer t_symbuf
        inline unsigned int get_bufsize() const {return bufsize; }
        //! maximum number of characters stored in an entry
        inline unsigned int get_max_entry_size() const {return max_entry_size; }
        /*! determine the actual current number of entries and amount of allocated memory,
            The statistics update is necessary because unused entries might have been
            deallocated after the last shared_ptr referencing them vanished.
            The whole memory estimate is anyway a) inaccurate because the operator sizeof
            does not count additional memory allocated from within the smart pointers,
            and b) ambiguous because it is unclear whether the handed out shared pointers
            should also be included in the memory accountance
            a) could be solved by overloading the operators new/delete with a custom allocator
               that exactly counts the allocated memory, but that would decrease performance. */
        void update_stats();
        //! the number of elements stored
        inline size_t get_nr()
        {update_stats(); return nr; }
        //! the number of elements stored
        inline size_t get_nbuf()
        {update_stats(); return nbuf; }
        /*! get the amount of memory currently occupied by the symbol table
            the contribution due to the string entries might not be totally accurate */
        inline size_t get_total_mem()
        {update_stats(); return mem; }
        //! check whether an entry of the symbol table exists
        bool contains_entry(const std::string &s) const
        {return contains_entry(s.c_str()); }
        bool contains_entry(const char * const s) const
        {
            std::lock_guard<std::mutex> lock(this->mtx);
            return (nonconst_get_entry(s) != nullptr);
        }
        /* ! condense information from get_topology, determine how often a given string
         *   is referenced by how many and distinct shared_ptr sets */
        void get_string_stats(const char * const s, size_t *ninst, size_t *nptr) const;
        /*! rename an entry and return true if the entry existed previously
            and false if it was created newly */
        bool check_rename_entry(const std::string &before, const std::string &after)
        {return check_rename_entry(before.c_str(), after.c_str()); }
        bool check_rename_entry(const char * const before, const char * const after)
        {return nonconst_rename_entry(before, after) != nullptr; }
    private:
        //!  only shared_string is granted direct access to the class members
        friend class shared_string;
        // private member functions
        /*! trim whitespace at the beginning of in and any chars beyond
            the maximum string length maxlen and return the trimmed string */
        std::string trim_string(const std::string &s, const size_t maxlen) const
        {return trim_string(s.c_str(), maxlen); }
        std::string trim_string(const char *s1, const size_t maxlen) const;
        //! create a new entry to be attached subsequently to a t_symbuf
        std::shared_ptr<std::shared_ptr<t_symbol> > construct_entry(const char * const in, const size_t maxlen, size_t *tmp_mem);
        //! private versions of get_entry returning a pointer to a non-const entry
        std::shared_ptr<std::shared_ptr<t_symbol> > nonconst_get_entry(const std::string &s) const {return nonconst_get_entry(s.c_str()); }
        std::shared_ptr<std::shared_ptr<t_symbol> > nonconst_get_entry(const char *s) const;
        //! private versions of insert_entry returning a pointer to a non-const entry
        std::shared_ptr<std::shared_ptr<t_symbol> > nonconst_insert_entry(const std::string &s){return nonconst_insert_entry(s.c_str()); }
        std::shared_ptr<std::shared_ptr<t_symbol> > nonconst_insert_entry(const char * const s);
        //! private versions of rename_entry returning a pointer to a non-const entry
        std::shared_ptr<std::shared_ptr<t_symbol> > nonconst_rename_entry(const std::string &before, const std::string &after){return nonconst_rename_entry(before.c_str(), after.c_str()); }
        std::shared_ptr<std::shared_ptr<t_symbol> > nonconst_rename_entry(const char * const before, const char * const after);
        /*! t_symtab is notified whenever a shared_string called its destructor
            when more shared_strings died than there are entries in a symbol buffer,
            clean() tries to reduce the occupied amount of memory by defragmenting
            buffers and deallocating empty buffers */
        void track_orphan_entries()
        {
            // FIXME: would a less stringent memory order suffice, e.g., memory_order_consume?
            if (n_orphan_entries.fetch_add((size_t)1, std::memory_order_seq_cst) > bufsize)
            {
                n_orphan_entries.store((size_t)0, std::memory_order_seq_cst);
                clean();
            }
        }
        //! determine pointer topology
        void get_topology(std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> > &topology) const;
        // private member variables
        //! entries per buffer class t_symbuf
        const unsigned int        bufsize;
        //! maximum size of an entry
        const unsigned int        max_entry_size;
        //! the number of symbol buffers in the linked list
        size_t                    nbuf;
        //! the number of symbol table entries
        size_t                    nr;
        //! pointer to the first block of entries enclosed in a t_symbuf structure
        std::shared_ptr<t_symbuf> symbuf;
        //! the memory occupied by our symtab
        std::atomic<size_t>       mem;
        //! the number of shared_str that called their destructor since the last clean()
        std::atomic<size_t>       n_orphan_entries;
        /*! indicates whether write access to the object is possible
            FIXME: does the use of std::mutex and std::thread create
            problems when other threading solutions are used even though we
            do not actually create threads ?*/
        mutable std::mutex mtx;
};

/// \endcond DEV
//! \}

}      // namespace gmx

#endif // symtab_cpp11_h
