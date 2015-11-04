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
    \inlibraryapi

   \brief
    Class Symtab stores std::strings non-redundantly (thread-safe) and
    trimmed to the chars actually used. Global replacement is possible, e.g.,
    for renaming all residues of a certain type in GROMACS.

    The class is intended as C++ replacement for the C implementation of the
    symbol table "Symtab" for convenient use just like a std::string and
    automatic, clean memory management.

    Only suitable for a limited number of unique entries/strings because the
    entries are stored unsorted, not allowing for time-saving lexical search.
    -> entry retrieval time grows linearly with the number of unique entries.

    The entries are not modifiable directly but only via dedicated member
    functions/operators to ensure internal consistency.

    The symbol table contains a linked list of symbol buffers Symbuf.
    Each Symbuf encloses a set of string entries defined by bufsize.
    Each string entry is saved as array of chars. Excessively long strings
    are truncated to a maximum length defined via maxEntrySize_.

    C++11 features are used for thread safety and improved memory management.
    C++11 smart pointers are used to keep track of the number of times
    entries are referenced and for automatic disposal of unreferenced entries.
    An entry ceases to exist as soon as no more external shared_ptr references
    are held.

   \author
    R. Thomas Ullmann <tullman@gwdg.de>

   \copyright
    GROMACS license

   \date Sep 2015
 */
#ifndef GMX_UTILITY_SYMTAB_CPP_H
#define GMX_UTILITY_SYMTAB_CPP_H

#include <cstring>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include <atomic>
#include <mutex>
#include <thread>


namespace gmx
{

/*!
   \class Symbol

   \brief
    Object that stores a single string.

   \ingroup module_utility
   \libinternal
 */
class Symbol
{
    public:
        //! the constructor
        Symbol()
        {
            ptr_.store(nullptr);
            group_.store(-1);
            cnt_.store(1);
        }
        //! the destructor
        ~Symbol()
        {
            if (ptr_.load(std::memory_order_seq_cst) != nullptr)
            {
                const char* tmp_p = ptr_.exchange(nullptr);
                if (tmp_p != nullptr)
                {
                    delete [] tmp_p;
                }
            }
        }

    private:
        //! only classes constituting shared_string can access the class members
        friend class Symbuf;
        friend class Symtab;
        friend class shared_string;

        // private member functions:

        //! determine how often the string/symbol is referenced, alias of getRef for consistency with shared_ptr.use_count()
        size_t use_count() const { return getRef(); }
        //! determine how often the string/symbol is referenced
        size_t getRef() const { return cnt_.load(std::memory_order_seq_cst); }
        //! set the reference count
        //! \param[in]   val   new reference count
        void   setRef(const size_t &val) { cnt_.store(val, std::memory_order_seq_cst); }
        //! increment the reference count
        size_t incrRef()
        { return (cnt_.fetch_add((size_t) 1, std::memory_order_seq_cst) + (size_t)1); }
        //! decrement the reference count
        size_t decrRef()
        { return (cnt_.fetch_add((size_t)-1, std::memory_order_seq_cst) - (size_t)1); }

        //! retrieve a pointer to the stored string
        char* getPtr() { return ptr_.load(std::memory_order_seq_cst); }
        //! retrieve a pointer to the stored string, alias of getPtr for consistency with share_ptr.get()
        char* get() { return getPtr(); }
        //! let the internal pointer reference a new string
        //! \param[in]   new_ptr   input pointer
        void setPtr(char* const new_ptr)
        { return ptr_.store(new_ptr, std::memory_order_seq_cst); }
        //! let the internal pointer reference a new string
        //! \param[in]   new_ptr   input pointer
        char* exchangePtr(char* const new_ptr)
        { return ptr_.exchange(new_ptr, std::memory_order_seq_cst); }

        //! retrieve the string group ID
        long int getGroup() { return group_.load(std::memory_order_seq_cst); }
        //! set a new string group ID
        //! \param[in]   new_group   input pointer
        void setGroup(const long int new_group)
        { return group_.store(new_group, std::memory_order_seq_cst); }
        //! exchange the string group ID with new_group and return the previous value
        //! \param[in]   new_group   input group ID
        int exchangeGroup(const long int new_group)
        { return group_.exchange(new_group, std::memory_order_seq_cst); }

        // private member variables:
        //! pointer to the actual string
        std::atomic<char*>  ptr_;
        //! optional string group, values < 0 mean no particular group membership,
        //! values >= 0 indicate a speparate string group, global replacement of
        //! strings can be restricted to a particular group in this way
        std::atomic<long int> group_;
        //! 2nd level reference count:
        //! number of distinct sets of shared_ptr "spp" referencing
        //! the shared pointer that holds the atomic string
        //! spp_{1,...,n} -> sp -> Symbol
        std::atomic<size_t> cnt_;
};

/*!
   \class Symbuf

   \brief
    A buffer object that stores/references a preset number of symbol objects.

   \ingroup module_utility
   \libinternal
 */
class Symbuf
{
    public:
        //! the default constructor
        Symbuf() : bufsize_(0), buf_(nullptr), prev_(), next_(nullptr)
        {
        }
        //! constructor with buffer preallocation for bs entries
        //! \param[in]   bs   buffer size
        Symbuf(const size_t bs) : bufsize_(bs), buf_(nullptr), prev_(), next_(nullptr)
        {
            allocSymbuf();
        }
        //! explicitly state that no copying is allowed
        Symbuf(const Symbuf &) = delete;
        //! the destructor
        ~Symbuf()
        {
            // by default, we leave the string entries intact for further use
            // the shared_ptr will take care of their deallocation later on
            freeSymbuf(false);
            next_ = nullptr;
        }
    private:
        //!  only Symtab is granted direct access to the class members
        friend class Symtab;

        // member variables

        //! maximum number of stored chars
        size_t bufsize_;
        //! pointer to the stored data (array of smart pointers to strings)
        //! use weak_ptr to ensure that a string is automatically destroyed
        //! as soon as there is no more external shared_ptr referencing it
        std::weak_ptr<std::shared_ptr<Symbol> > *buf_;
        //! pointer to the previous Symbuf, weak_ptr to prevent cyclic blockage of deallocation
        std::weak_ptr<Symbuf>                    prev_;
        //! pointer to the next Symbuf
        std::shared_ptr<Symbuf>                  next_;

        // member functions

        //! compare two strings
        //! \param[in]   s1   first  string to compare
        //! \param[in]   s2   second string to compare
        bool matchStrings(const std::string &s1, const std::string &s2) const {return matchStrings(s1.c_str(), s2.c_str()); }
        //! compare two strings
        //! \param[in]   s1   first  string to compare
        //! \param[in]   s2   second string to compare
        bool matchStrings(const char * const s1, const char * const s2) const;

        //! dellocate the Symbuf member array, but leave the referenced strings intact
        void freeSymbuf(){ return freeSymbuf(false); }
        //! deallocate all stored string entries
        //! \param[in]   dealloc_entries   true: deallocate all string entries, false: keep all string entries intact
        void freeSymbuf(const bool &dealloc_entries)
        {
            if (buf_ != nullptr)
            {
                // on request, enforce deletion of the stored strings
                // even if there are still shared_ptrs referencing them
                if (dealloc_entries)
                {
                    for (size_t i = 0; i < bufsize_; ++i)
                    {
                        const char* tmp_cptr = (buf_[i].lock()->get())->exchangePtr(nullptr);
                        delete [] tmp_cptr;
                    }
                }
                for (size_t i = 0; i < bufsize_; ++i)
                {
                    // detach the weak_ptr from its target
                    buf_[i].reset();
                }
                ::operator delete (buf_);
                buf_ = nullptr;
            }
            bufsize_ = 0;
        }
        //! allocate memory for the pointers that will reference the string entries
        //! returns the amount of memory occupied by *this buffer after allocation
        size_t allocSymbuf()
        {
            try // allocate memory for the weak_ptrs
            {
                const size_t this_array_size = sizeof(std::weak_ptr<std::shared_ptr<Symbol> >) * bufsize_;
                buf_ = static_cast< std::weak_ptr<std::shared_ptr<Symbol> >* > (::operator new(this_array_size));
                for (size_t i = 0; i < bufsize_; ++i)
                {
                    // create the weak_ptrs in the allocated memory region
                    ::new (&buf_[i]) std::weak_ptr<std::shared_ptr<std::shared_ptr<Symbol> > >();
                    if (!buf_[i].expired())
                    {
                        buf_[i].reset();
                    }
                }
                size_t mem = sizeof(Symbuf) + this_array_size;
                return mem;
            }
            catch (std::bad_alloc &)
            {
                fprintf (stderr, "Error in %s::allocSymbuf():",  typeid(*this).name());
                fprintf (stderr, "Insufficient memory for allocating new symbol table buffer!\n");
                throw;
            }
            return 0;
        }

        //! retrieve a string entry
        //! \param[in]   s       query string
        //! \param[in]   group   group id
        std::shared_ptr<std::shared_ptr<Symbol> > getEntry(const char * const s, const long int group) const;

        //! retrieve a string entry
        //! \param[in]   s       query string
        std::shared_ptr<std::shared_ptr<Symbol> > getEntry(const char * const s) const
        { return getEntry(s, (long int)-1); }

        //! check existence of a string entry
        //! \param[in]   s       query string
        //! \param[in]   group   group id
        bool containsEntry(const char * const s, const long int group) const
        { return (getEntry(s, group) != nullptr); }

        //! check existence of a string entry
        //! \param[in]   s       query string
        bool containsEntry(const char * const s) const
        { return (getEntry(s) != nullptr); }

        //! store a string entry
        //! \param[in]   in   input Symbol referenced by two layers of shared_ptrs to enable global renaming
        bool insertEntry(std::shared_ptr<std::shared_ptr<Symbol> > in);

        //! get number of strings stored and an estimate of the total amount of memory occupied by this Symbuf
        //! \param[out]   nr    number of strings
        //! \param[out]   mem   estimate of the total memory usage in byte
        void getStats(unsigned int * const nr, size_t * const mem);
        //! determine the pointer topology of this Symbuf
        //! \param[out]   topology   first hierarchy level: string groups, second level: string, third level:
        //!                          map of shared pointers from distinct shared_ptr sets referencing the string,
        //!                          number of active shared_ptr instances of the respective set
        void getTopology(std::map<long int, std::map<std::string, std::map<std::shared_ptr<Symbol>*, size_t> > > &topology) const;
};

/*!
   \class Symtab

   \brief
    non-redundant storage of strings with std::string-like usage,
    enables global replacement

   \throws std::bad_alloc if out of memory

   \throws std::logic_error in case of an (obviously) broken implementation

   A C++ replacement for the C implementation of the symbol table "Symtab".
   Uses a linked list of storage objects of class Symbuf, which reference multiple
   Symbol objects. Each individual string is stored by one of these Symbol objects.

   \ingroup module_utility
   \inlibraryapi
 */
class Symtab
{
    public:
        //! the default constructor
        Symtab() : bufsize_(81), maxEntrySize_(1024), nbuf_(0),
                   nr_(0), symbuf_(nullptr), mem_(sizeof(*this))
        {
        }
        //! a constructor that allows custom values for buffer size and maximum entry length
        //! \param[in]  bs   buffer size (number of string entries that can be stored)
        //! \param[in]  mes  maximum string length of an entry above which the string is truncated
        Symtab(const size_t &bs, const size_t &mes) : bufsize_(bs), maxEntrySize_(mes), nbuf_(0),
                                                      nr_(0), symbuf_(nullptr), mem_(sizeof(*this))
        {
        }
        //! explicitly state that no copying is allowed
        Symtab(const Symtab &) = delete;

        //! the destructor
        ~Symtab()
        {
            // should all happen by itself thanks to the shared_pointer linkage of the Symbufs
            // freeSymtab();
            // symbuf_.reset();
        }
        //! by default, leave the strings intact, shared_ptrs will free them on destruction
        void freeSymtab(){ return freeSymtab(false); }

        //! deallocate the symbol table and if desired also the stored string entries
        //! \param[in]   dealloc_entries   false: strings managed by the handed out shared_ptrs stay intact,
        //!                                true:  strings are deallocated
        void freeSymtab (const bool &dealloc_entries)
        {
            std::shared_ptr<Symbuf> tmp_buf(symbuf_);
            // go along the linked list of symbufs deleting all entries,
            // or only the weak references to them if !dealloc_entries
            while (tmp_buf != nullptr)
            {
                tmp_buf->freeSymbuf(dealloc_entries);
                tmp_buf = tmp_buf->next_;
            }
            // initiate the successive calling of dtors along the linked list of Symbufs
            symbuf_.reset();
        }

        //! \brief defragment Symbufs by moving entries towards the beginning of
        //!        the linked list and remove left over empty Symbuf objects
        void clean();

        //! maximum number of elements stored per symbol buffer Symbuf
        inline unsigned int getBufsize() const { return bufsize_; }

        //! maximum number of characters stored in an entry
        inline unsigned int getMaxEntrySize() const { return maxEntrySize_; }

        /*! \brief determine the actual current number of entries and the approximate total
                   amount of allocated memory.

            The statistics update is necessary because unused entries might have been
            deallocated after the last shared_ptr referencing them vanished.
            The whole memory estimate is anyway a) inaccurate because the operator sizeof
            does not count additional memory allocated from within the smart pointers,
            and b) ambiguous because it is unclear whether the handed out shared pointers
            should also be included in the memory accountance
            a) could be solved by overloading the operators new/delete with a custom allocator
               that exactly counts the allocated memory, but that would decrease performance
               and is not particularly popular practice */
        void updateStats();

        //! the number of elements stored
        inline size_t getNr()
        { updateStats(); return nr_; }

        //! the number of elements stored
        inline size_t getNbuf()
        { updateStats(); return nbuf_; }

        //! get the amount of memory currently occupied by the symbol table
        //!  the contribution due to the string entries might not be totally accurate
        inline size_t getTotalMem()
        { updateStats(); return mem_; }

        //! check whether an entry of the symbol table exists
        //! \param[in]   s   query string
        bool containsEntry(const std::string &s) const
        { return containsEntry(s.c_str()); }
        //! determine whether an entry exists
        //! \param[in]   s   query string
        bool containsEntry(const char * const s) const
        {
            return (containsEntry(s, (long int)-1));
        }

        //! check whether an entry of the symbol table exists
        //! \param[in]   s       query string
        //! \param[in]   group   string group ID
        bool containsEntry(const std::string &s, const long int group) const
        { return containsEntry(s.c_str(), group); }
        //! determine whether an entry exists
        //! \param[in]   s       query string
        //! \param[in]   group   string group ID
        bool containsEntry(const char * const s, const long int group) const
        {
            std::lock_guard<std::mutex> lock(mtx_);
            return (nonconstGetEntry(s, group) != nullptr);
        }

        //! rename an entry, return true on success or false on failure
        //! \param[in]   before   string entry to be replaced
        //! \param[in]   after    replacement string
        bool checkRenameEntry(const std::string &before, const std::string &after)
        { return checkRenameEntry(before.c_str(), after.c_str()); }

        //! rename an entry, return true on success or false on failure
        //! \param[in]   before   string entry to be replaced
        //! \param[in]   after    replacement string
        bool checkRenameEntry(const char * const before, const char * const after)
        { return checkRenameEntry(before, (long int)-1, after, (long int)-1); }

        //! rename an entry, return true on success or false on failure
        //! \param[in]   before         string entry to be replaced
        //! \param[in]   group_before   current string group ID
        //! \param[in]   after          replacement string
        //! \param[in]   group_after    new string group ID
        bool checkRenameEntry(const std::string &before, const long int group_before,
                              const std::string &after,  const long int group_after)
        { return checkRenameEntry(before.c_str(), group_before, after.c_str(), group_after); }

        //! rename an entry, return true on success or false on failure
        //! \param[in]   before         string entry to be replaced
        //! \param[in]   group_before   current string group ID
        //! \param[in]   after          replacement string
        //! \param[in]   group_after    new string group ID
        bool checkRenameEntry(const char * const before, const long int group_before,
                              const char * const after,  const long int group_after)
        { return nonconstRenameEntry(before, group_before, after, group_after) != nullptr; }

        //! condense information from getTopology, determine how often a given string
        //! is referenced by how many and distinct shared_ptr sets
        //! \param[in]    s       query string
        //! \param[out]   ninst   number of shared_string instances using the query string
        //! \param[out]   nptr    number of distinct outer shared_ptr sets using the query
        //!                       string, usually 1, but renaming can lead to multiple sets
        //!                       of outer pointers referencing the string
        void getStringStats(const char * const s, size_t *ninst, size_t *nptr) const
        { return getStringStats(s, (long int)-1, ninst, nptr); }
        //! condense information from getTopology, determine how often a given string
        //! is referenced by how many and distinct shared_ptr sets
        //! \param[in]    s       query string
        //! \param[in]    group   string group of the query string
        //! \param[out]   ninst   number of shared_string instances using the query string
        //! \param[out]   nptr    number of distinct outer shared_ptr sets using the query
        //!                       string, usually 1, but renaming can lead to multiple sets
        //!                       of outer pointers referencing the string
        void getStringStats(const char * const s, long int group, size_t *ninst, size_t *nptr) const;

    private:
        //!  only shared_string is granted direct access to the class members
        friend class shared_string;

        // private member functions

        //! trim whitespace at the beginning of in and any chars beyond
        //! the maximum string length maxlen and return the trimmed string
        //! \param[in]   s        input string
        //! \param[in]   maxlen   maximum string length
        std::string trimString(const std::string &s, const size_t maxlen) const
        { return trimString(s.c_str(), maxlen); }
        //! trim whitespace at the beginning of in and any chars beyond
        //! the maximum string length maxlen and return the trimmed string
        //! \param[in]   s        input string
        //! \param[in]   maxlen   maximum string length
        std::string trimString(const char *s, const size_t maxlen) const;

        //! create a new entry to be attached subsequently to a Symbuf
        //! \param[in]       in         string to be stored
        //! \param[in]       group_in   string group to which the entry will be added
        //! \param[in]       maxlen     maximum string length stored, rest trimmed
        //! \param[in,out]   tmp_mem    memory occupied by the Symtab, will be updated
        std::shared_ptr<std::shared_ptr<Symbol> > constructEntry(const char * const in, const long int group_in, const size_t maxlen, size_t *tmp_mem);

        //! private version of getEntry returning a pointer to a non-const entry
        //! if the query string is found and nullptr otherwise
        //! \param[in]       s         string to be retrieved
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstGetEntry(const std::string &s) const
        { return nonconstGetEntry(s.c_str()); }
        //! private version of getEntry returning a pointer to a non-const entry
        //! if the query string is found and nullptr otherwise
        //! \param[in]       s         string to be retrieved
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstGetEntry(const char *s) const
        { return nonconstGetEntry(s, (long int)-1); }
        //! private version of getEntry returning a pointer to a non-const entry
        //! if the query string is found and nullptr otherwise
        //! \param[in]       s         string to be retrieved
        //! \param[in]       group     string group in which to search for s
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstGetEntry(const std::string &s, const long int group) const
        { return nonconstGetEntry(s.c_str(), group); }
        //! private version of getEntry returning a pointer to a non-const entry
        //! if the query string is found and nullptr otherwise
        //! \param[in]       s         string to be retrieved
        //! \param[in]       group     string group in which to search for s
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstGetEntry(const char *s, const long int group) const;

        //! private version of insertEntry returning a pointer to a non-const entry
        //! \param[in]       s         string to be retrieved
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstInsertEntry(const std::string &s)
        { return nonconstInsertEntry(s.c_str()); }
        //! private version of insertEntry returning a pointer to a non-const entry
        //! \param[in]       s         string to be retrieved
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstInsertEntry(const char * const s)
        { return nonconstInsertEntry(s, (long int)-1); }
        //! private version of insertEntry returning a pointer to a non-const entry
        //! \param[in]       s         string to be retrieved
        //! \param[in]       group     string group to add s to
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstInsertEntry(const std::string &s, const long int group)
        { return nonconstInsertEntry(s.c_str(), group); }
        //! private version of insertEntry returning a pointer to a non-const entry
        //! \param[in]       s         string to be retrieved
        //! \param[in]       group     string group to add s to
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstInsertEntry(const char * const s, const long int group);

        //! private version of rename_entry returning a pointer to a non-const entry
        //! \param[in]   before   string entry to be replaced
        //! \param[in]   after    replacement string
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstRenameEntry(const std::string &before, const std::string &after)
        { return nonconstRenameEntry(before.c_str(), after.c_str()); }
        //! private version of rename_entry returning a pointer to a non-const entry
        //! \param[in]   before   string entry to be replaced
        //! \param[in]   after    replacement string
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstRenameEntry(const char * const before, const char * const after)
        { return nonconstRenameEntry(before, (long int)-1, after, (long int)-1); }
        //! private version of rename_entry returning a pointer to a non-const entry
        //! \param[in]   before         string entry to be replaced
        //! \param[in]   group_before   current string group ID
        //! \param[in]   after          replacement string
        //! \param[in]   group_after    new string group ID
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstRenameEntry(std::string &before, const long int group_before,
                                                                      std::string &after,  const long int group_after)
        { return nonconstRenameEntry(before.c_str(), group_before, after.c_str(), group_after); }
        //! private version of rename_entry returning a pointer to a non-const entry
        //! \param[in]   before         string entry to be replaced
        //! \param[in]   group_before   current string group ID
        //! \param[in]   after          replacement string
        //! \param[in]   group_after    new string group ID
        std::shared_ptr<std::shared_ptr<Symbol> > nonconstRenameEntry(const char * const before, const long int group_before,
                                                                      const char * const after,  const long int group_after);

        //! determine pointer topology for debugging purposes
        //! \param[out]   topology   first hierarchy level: string groups, second level: string, third level:
        //!                          map of shared pointers from distinct shared_ptr sets referencing the string,
        //!                          number of active shared_ptr instances of the respective set
        void getTopology(std::map<long int, std::map<std::string, std::map<std::shared_ptr<Symbol>*, size_t> > > &topology) const;

        // private member variables

        //! entries per buffer class Symbuf
        const unsigned int        bufsize_;
        //! maximum size of an entry
        const unsigned int        maxEntrySize_;
        //! the number of symbol buffers in the linked list
        size_t                    nbuf_;
        //! the number of symbol table entries
        size_t                    nr_;
        //! pointer to the first block of entries enclosed in a Symbuf structure
        //! additional Symbuf objects may follow in a linked list
        std::shared_ptr<Symbuf>   symbuf_;
        //! the memory occupied by our symtab
        std::atomic<size_t>       mem_;
        //! indicates whether write access to the object is possible
        mutable std::mutex        mtx_;
};


}      // namespace gmx

#endif // end header
