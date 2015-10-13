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

   The entries of the symbol table are not modifiable directly but only via
   dedicated member functions/operators to ensure internal consistency.

   The symbol table contains a linked list of symbol buffers t_symbuf.
   Each t_symbuf encloses a set of string entries defined by bufsize.
   Each string entry is saved as array of chars. Excessively long strings
   are truncated to a maximum length defined via max_entry_size.

   Mutex and atomics from thread-MPI are used for thread safety.

   The boost smart pointers included with GROMACS are used to keep track of
   the number of times entries are referenced and for automatic disposal of
   unreferenced entries.

   \author
    R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \copyright
    GROMACS license

   \date Feb 2015
 */

#ifndef gmx_symtab_cpp98_h
#define gmx_symtab_cpp98_h

#include <string>
#include <cstring>
#include <map>
#include <stdexcept>

// change as of commit d56c6f64ef2ba1ec31ec495c5204adbad24ab9ee
// still gmx::Mutex and gmx_lock_guard wrap the ThreadMPI
// implementations, but may be exchanged behind the scenes
//#include "external/thread_mpi/include/thread_mpi/mutex.h"
#include "gromacs/utility/mutex.h"
#include "external/thread_mpi/include/thread_mpi/atomic.h"
#include "external/boost/boost/shared_ptr.hpp"

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

/*!
   \class t_symbol shared_string.h "gromacs/data_structures/shared_string.h"

   \brief
    storage of a "symbol", Here, symbol refers to a std::string of characters.

   \param[in] char* pointer to a char std::string to be stored or looked up.

   \returns shared_ptr<shared_ptr<t_symbol> > upon inserting or retrieving an entry,
            where t_symbol stores the std::string

   \throws std::bad_alloc if out of memory

   \throws std::logic_error in case of an (obviously) broken implementation
 */
class t_symbol
{
    public:
        t_symbol()
        {
            tMPI_Atomic_ptr_set(&ptr, NULL);
            tMPI_Atomic_set(&cnt, (int)1);
        }
        ~t_symbol()
        {
            void *tp = tMPI_Atomic_ptr_get(&ptr);
            if (tp != NULL)
            {
                const char* tmp_p = static_cast<char*>(tp);
                // there is no atomic pointer swap function without comparison in tMPI
                if (tmp_p != NULL && tMPI_Atomic_ptr_cas(&ptr, tp, NULL) != 0)
                {
                    delete [] tmp_p;
                }
                else
                {
                    fprintf(stderr, "Error in destructor of %s:\n", typeid(*this).name());
                    fprintf(stderr, "Pointer value changed while trying to deallocate the stored symbol.\n");
                    throw std::logic_error("Internal error in t_symtab.");
                }
            }
        }
    private:
        //! only classes constituting shared_string can access the class members
        friend class t_symbuf;
        friend class t_symtab;
        friend class shared_string;
        friend t_symbol* symtab_inner_allocator(const size_t &len);
        friend boost::shared_ptr<t_symbol>* symtab_outer_allocator(const size_t &len);
        friend void symtab_inner_deleter(t_symbol* &p);
        friend void symtab_outer_deleter(boost::shared_ptr<t_symbol>* &p);
        // private member functions:
        size_t use_count() const
        {
            return get_ref();
        }
        size_t get_ref() const
        {
            return static_cast<size_t>(tMPI_Atomic_get(&cnt));
        }
        void   set_ref(const size_t &val)
        {
            tMPI_Atomic_set(&cnt, val);
        }
        size_t incr_ref()
        {
            return static_cast<size_t>(tMPI_Atomic_add_return(&cnt, (int) 1));
        }
        size_t decr_ref()
        {
            return static_cast<size_t>(tMPI_Atomic_add_return(&cnt, (int)-1));
        }
        char* get_ptr()
        {
            return static_cast<char*>(tMPI_Atomic_ptr_get(&ptr));
        }
        void set_ptr(char * new_ptr)
        {
            tMPI_Atomic_ptr_set(&ptr, static_cast<void*>(new_ptr));
        }
        char* exchange_ptr(char * new_ptr)
        {
            void *tp;
            do
            {
                tp = tMPI_Atomic_ptr_get(&ptr);
            }
            while (tMPI_Atomic_ptr_cas(&ptr, tp, new_ptr) == 0);

            return static_cast<char*>(tp);
        }
        // private member variables:
        //! pointer to the actual std::string
        tMPI_Atomic_ptr_t ptr;
        // 2nd level reference count:
        // number of distinct sets of shared_ptr "spp" referencing
        // the shared pointer that holds the atomic std::string
        // spp_{1,...,n} -> sp -> t_symbol
        //! 2nd level reference count
        tMPI_Atomic_t cnt;
};

/*!
   \class t_symbuf shared_string.h "gromacs/data_structures/shared_string.h"

   \brief
    storage of t_symbol objects for class t_symtab

   \param[in] char* pointer to a char std::string to be stored or looked up.

   \returns shared_ptr<shared_ptr<t_symbol> > upon inserting or retrieving an entry,
            where t_symbol stores the std::string

   \throws std::bad_alloc if out of memory
 */
class t_symbuf
{
    public:
        //! the default constructor
        t_symbuf() : prev()
        {
            bufsize   = 0;
            buf       = NULL;
            next.reset();
        }
        //! the destructor
        ~t_symbuf()
        {
            // by default, we leave the std::string entries intact for further use
            // the shared_ptr will take care of their deallocation later on
            free_symbuf(false);
            next.reset();
        }
    private:
        //!  only t_symtab is granted direct access to the class members
        friend class t_symtab;
        friend t_symbol* symtab_inner_allocator(const size_t &len);
        friend boost::shared_ptr<t_symbol>* symtab_outer_allocator(const size_t &len);
        friend void symtab_inner_deleter(t_symbol* &p);
        friend void symtab_outer_deleter(boost::shared_ptr<t_symbol>* &p);
        // member variables
        //! maximum number of stored chars
        size_t bufsize;
        /*! \brief pointer to the stored data (array of smart pointers to std::strings)

            In absence of weak_ptr, make sure that the deallocator takes care that
            a std::string is automatically destroyed as soon as there is no more external
            shared_ptr referencing it */
        boost::shared_ptr<boost::shared_ptr<t_symbol> > *buf;
        //! pointer to the next t_symbuf
        boost::shared_ptr<t_symbuf>                      next;
        /*! \brief pointer to the previous t_symbuf, In absence of weak_ptr,
            the destructor has to avoid cyclic blockage of deallocation. */
        boost::shared_ptr<t_symbuf>                     prev;
        // member functions
        std::string trim_string(const char * const s, const size_t maxlen) const;
        bool match_strings(const std::string &s1, const std::string &s2) const {return match_strings(s1.c_str(), s2.c_str()); }
        bool match_strings(const char * const s1, const char * const s2) const;
        void free_symbuf(){ return free_symbuf(false); }
        void free_symbuf(const bool &dealloc_entries)
        {
            if (buf != NULL)
            {
                // on request, enforce deletion of the stored std::strings
                // even if there are still shared_ptrs referencing them
                if (dealloc_entries)
                {
                    for (size_t i = 0; i < bufsize; ++i)
                    {
                        const char* tmp_cptr = ((buf[i])->get())->exchange_ptr(NULL);
                        delete [] tmp_cptr;
                    }
                }
                for (size_t i = 0; i < bufsize; ++i)
                {
                    buf[i].reset();
                }
                ::operator delete (buf);
                buf = NULL;
            }
            bufsize = 0;
        }
        void remove_orphan_entries()
        {
            if (buf != NULL)
            {
                for (size_t i = 0; i < bufsize; ++i)
                {
                    if (buf[i].unique())
                    {
                        buf[i].reset();
                    }
                }
            }
        }
        void alloc_symbuf(size_t * const mem)
        {
            try // allocate memory for the shared_ptrs
            {
                const size_t this_array_size = sizeof(boost::shared_ptr<boost::shared_ptr<t_symbol> >) * bufsize;
                buf = static_cast<boost::shared_ptr<boost::shared_ptr<t_symbol> >* > (::operator new(this_array_size));
                for (size_t i = 0; i < bufsize; ++i)
                {
                    // create the shared_ptrs in the allocated memory region
                    ::new (&buf[i]) boost::shared_ptr<boost::shared_ptr<boost::shared_ptr<t_symbol> > >();
                    buf[i].reset();
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
        // check existence of a std::string entry
        bool contains_entry(const char * const s) const { return (get_entry(s) != NULL); }
        //! retrieve a std::string entry
        boost::shared_ptr<boost::shared_ptr<t_symbol> > get_entry(const char * const s) const;
        //! store a std::string entry
        bool insert_entry(boost::shared_ptr<boost::shared_ptr<t_symbol> > in);
        //! get number of std::strings stored and an estimate of the total amount of memory occupied by this t_symbuf
        void get_stats(size_t * const nr, size_t * const mem);
        void get_topology(std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> > &topology) const;
};

/*!
   \class t_symtab shared_string.h "gromacs/data_structures/shared_string.h"

   \brief
    non-redundant storage of std::strings for class shared_string,
    enables global replacement.
    The t_symbol object holds a pointer to the actual std::string.
    Each t_symbuf holds a number (bufsize) of pointers to t_symbol objects.

   \param[in] char* pointer to a char std::string to be stored or looked up.

   \param[in] std::string a std::string to be stored or looked up.

   \returns shared_ptr<shared_ptr<t_symbol> > upon inserting, retrieving, modifying an entry,
            where t_symbol stores the std::string

   \throws std::bad_alloc if out of memory

   \throws std::logic_error in case of an (obviously) broken implementation
 */
class t_symtab
{
    public:
        //! the default constructor
        t_symtab() : bufsize(81), max_entry_size(1024), symbuf(), nbuf(0), nr(0), mem(sizeof(*this))
        {
            tMPI_Atomic_set(&n_orphan_entries, 0);
        }
        //! a constructor that allows custom values for buffer size and maximum entry length
        t_symtab(const unsigned int &bs, const unsigned int &mes) : bufsize(bs), max_entry_size(mes), symbuf(),
                                                                    nbuf(0), nr(0), mem(sizeof(*this))
        {
            tMPI_Atomic_set(&n_orphan_entries, 0);
        }
#if  __cplusplus >= 201103L
        //! explicitly state that no copying is allowed
        t_symtab(const t_symtab &) = delete;
#endif
        //! the destructor
        ~t_symtab()
        {
            //! with pointer prev in t_symbuf as weak_pointer this, would all happen by itself
            //! thanks to the shared_pointer linkage of the t_symbufs
            free_symtab();
        }
        //! by default, leave the std::strings intact, shared_ptrs will free them on destruction
        void free_symtab(){ return free_symtab(false); }
        // (dealloc_entries == false) -> std::strings managed by the handed
        // out shared_ptrs stay intact
        void free_symtab (const bool &dealloc_entries)
        {
            boost::shared_ptr<t_symbuf> tmp_buf(symbuf);
            // go along the linked list of symbufs deleting all entries,
            // or only the symbuf references to them if !dealloc_entries
            while (tmp_buf != NULL)
            {
                tmp_buf->free_symbuf(dealloc_entries);
                tmp_buf = tmp_buf->next;
                //! with pointer prev in t_symbuf as weak_pointer this would not be needed
                if (tmp_buf != NULL)
                {
                    tmp_buf->prev.reset();
                }
            }
            // initiate the successive calling of dtors along the linked list of t_symbufs
            symbuf.reset();
        }
        /*! defragment t_symbufs by moving entries towards the beginning of the linked list.
            remove left over empty t_symbuf structures */
        void clean();
        //! maximum number of elements stored per symbol buffer t_symbuf
        inline unsigned int get_bufsize() const { return bufsize; }
        //! maximum number of characters stored in an entry
        inline unsigned int get_max_entry_size() const { return max_entry_size; }
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
        { update_stats(); return nr; }
        //! the number of elements stored
        inline size_t get_nbuf()
        { update_stats(); return nbuf; }
        /*! get the amount of memory currently occupied by the symbol table
            the contribution due to the std::string entries might not be totally accurate */
        inline size_t get_total_mem()
        { update_stats(); return mem; }
        //! check whether an entry of the symbol table exists
        bool contains_entry(const std::string &s) const
        { return contains_entry(s.c_str()); }
        bool contains_entry(const char * const s) const
        {
            gmx::lock_guard<gmx::Mutex> lock(this->mtx);
            return (nonconst_get_entry(s) != NULL);
        }
        /* ! \brief condense information from get_topology, determine how often a given std::string
         *   is referenced by how many distinct shared_ptr sets */
        void get_string_stats(const char * const s, size_t *ninst, size_t *nptr) const;
        /*! \brief rename an entry and return true if the entry
            existed previously and false if it was created newly */
        bool check_rename_entry(const std::string &before, const std::string &after)
        { return check_rename_entry(before.c_str(), after.c_str()); }
        bool check_rename_entry(const char * const before, const char * const after)
        { return nonconst_rename_entry(before, after) != NULL; }
    private:
        //!  only shared_string is granted direct access to the class members
        friend class shared_string;
        friend t_symbol* symtab_inner_allocator(const size_t &len);
        friend boost::shared_ptr<t_symbol>* symtab_outer_allocator(const size_t &len);
        friend void symtab_inner_deleter(t_symbol* &p);
        friend void symtab_outer_deleter(boost::shared_ptr<t_symbol>* &p);
        // private member functions
        /*! trim whitespace at the beginning of in and any chars beyond
            the maximum std::string length maxlen and return the trimmed std::string */
        std::string trim_string(const std::string &s, const size_t maxlen) const
        {
            return trim_string(s.c_str(), maxlen);
        }
        std::string trim_string(const char *s1, const size_t maxlen) const;
        //! create a new entry to be attached subsequently to a t_symbuf
        boost::shared_ptr<boost::shared_ptr<t_symbol> > construct_entry(const char * const in, const size_t maxlen, size_t *tmp_mem);
        //! private versions of get_entry returning a pointer to a non-const entry
        boost::shared_ptr<boost::shared_ptr<t_symbol> > nonconst_get_entry(const std::string &s) const { return nonconst_get_entry(s.c_str()); }
        boost::shared_ptr<boost::shared_ptr<t_symbol> > nonconst_get_entry(const char *s) const;
        //! private versions of insert_entry returning a pointer to a non-const entry
        boost::shared_ptr<boost::shared_ptr<t_symbol> > nonconst_insert_entry(const std::string &s){ return nonconst_insert_entry(s.c_str()); }
        boost::shared_ptr<boost::shared_ptr<t_symbol> > nonconst_insert_entry(const char * const s);
        //! private versions of rename_entry returning a pointer to a non-const entry
        boost::shared_ptr<boost::shared_ptr<t_symbol> > nonconst_rename_entry(const std::string &before, const std::string &after){ return nonconst_rename_entry(before.c_str(), after.c_str()); }
        boost::shared_ptr<boost::shared_ptr<t_symbol> > nonconst_rename_entry(const char * const before, const char * const after);
        /*! t_symtab is notified whenever a shared_string called its destructor
            when more shared_strings died than there are entries in a symbol buffer,
            clean() tries to reduce the occupied amount of memory by defragmenting
            buffers and deallocating empty buffers */
        void track_orphan_entries()
        {
            //! FIXME(?) do these fences what I think they do? check! -- getting the whole
            //! shared_string assembly lock-free seems difficult and maybe not be worthwhile
            tMPI_Atomic_memory_barrier_acq();
            if (static_cast<size_t>(tMPI_Atomic_add_return(&n_orphan_entries, (int)1)) > bufsize)
            {
                tMPI_Atomic_set(&n_orphan_entries, (int)-1);
                clean();
                if (tMPI_Atomic_fetch_add(&n_orphan_entries, 0) != -1)
                {
                    throw std::logic_error("another thread interfered with t_symbuf::track_orphan_entries()");
                }
            }
            tMPI_Atomic_memory_barrier_rel();
        }
        //! determine pointer topology
        void get_topology(std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> > &topology) const;
        // private member variables
        //! entries per buffer class t_symbuf
        const unsigned int          bufsize;
        //! maximum size of an entry
        const unsigned int          max_entry_size;
        //! the number of shared_str that called their destructor since the last clean()
        tMPI_Atomic_t               n_orphan_entries;
        //! pointer to the first block of entries enclosed in a t_symbuf structure
        boost::shared_ptr<t_symbuf> symbuf;
        //! the number of symbol buffers in the linked list
        size_t                      nbuf;
        //! the number of symbol table entries
        size_t                      nr;
        //! the memory occupied by our symtab
        size_t                      mem;
        //! indicates whether write access to the object is possible
        mutable gmx::Mutex          mtx;
};

/// \endcond DEV

//! \}

}      // namespace gmx

#endif // gmx_symtab_cpp98_h
