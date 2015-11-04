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
    Implementation file needed for the symbol table Symtab

   \author
    R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \copyright
    GROMACS license

   \date Sep 2015
 */

#include "gmxpre.h"

#include "symtab_cpp.h"

namespace gmx
{

//=====================================================================================
//
// implementations of member functions for the symbol table objects
//
//=====================================================================================

//! function compares two strings, not considering wildcards, and
//! returns true if the strings are equal and false otherwise
//! leading and trailing whitespace is skipped
bool Symbuf::matchStrings(const char * const s1, const char * const s2) const
{
    if (s1 == nullptr || s2 == nullptr)
    {
        return false;
    }
    else if (*s1 == '\0' && *s2 == '\0')
    {
        return true;
    }
    const char *start1 = s1;
    const char *start2 = s2;

    // skip leading whitespace
    while (start1 != nullptr && (*start1 == ' ' || *start1 == '\t'))
    {
        start1++;
    }
    while (start2 != nullptr && (*start2 == ' ' || *start2 == '\t'))
    {
        start2++;
    }

    // find string ends
    const char *end1 = start1;
    const char *end2 = start2;
    // search the end of the string marked by the terminating '\0'
    while (end1 != nullptr && *end1 != '\0')
    {
        end1++;
    }
    while (end2 != nullptr && *end2 != '\0')
    {
        end2++;
    }
    // omit trailing whitespace, set ends to the last non-whitespace, non-terminating chars
    while (end1 != start1 && (*end1 == ' ' || *end1 == '\t' || *end1 == '\0'))
    {
        end1--;
    }
    while (end2 != start2 && (*end2 == ' ' || *end2 == '\t' || *end2 == '\0'))
    {
        end2--;
    }

    // compare all character positions
    while (start1 != end1 && start2 != end2)
    {
        if (*start1 != *start2)
        {
            return false;
        }
        start1++;
        start2++;
    }
    // check whether all characters could be compared
    if (*start1 != *start2 || start1 != end1 || start2 != end2)
    {
        return false;
    }

    return true;
}

bool Symbuf::insertEntry(std::shared_ptr<std::shared_ptr<Symbol> > in)
{
    // Search for a free pointer in the array buf.
    // If found, allocate memory for the new string entry
    // and store it in the newly allocated space
    for (size_t i = 0; i < bufsize_; ++i)
    {
        if (buf_[i].expired())
        {
            buf_[i] = in;
            return true;
        }
    }

    return false;
}

std::shared_ptr<std::shared_ptr<Symbol> > Symbuf::getEntry(const char * const s, const long int group) const
{
    if (s == nullptr)
    {
        return nullptr;
    }

    for (size_t i = 0; i < bufsize_; ++i)
    {
        if (buf_[i].expired())
        {
            continue;
        }
        else
        {
            std::shared_ptr<std::shared_ptr<Symbol> > tmp_ptr(buf_[i]);
            // negative group id means no particular group membership
            long int gid = group;
            if (group  < 0)
            {
                // negative group ID means no special group membership
                gid = -1;
            }
            if (matchStrings((**tmp_ptr).getPtr(), s) && (**tmp_ptr).getGroup() == gid)
            {
                return tmp_ptr;
            }
        }
    }

    return nullptr;
}

void Symbuf::getStats(unsigned int * const nr, size_t * const mem)
{
    (*nr)  = 0;
    (*mem) = sizeof(Symbuf) + sizeof(std::weak_ptr<std::shared_ptr<Symbol> >) * bufsize_;
    for (size_t i = 0; i < bufsize_; ++i)
    {
        if (!buf_[i].expired())
        {
            (*nr)++;
            size_t      len = 1 + std::strlen((**std::shared_ptr<std::shared_ptr<Symbol> >(buf_[i])).getPtr());
            long double one_over_secondary_ref_count = 1 /
                static_cast<long double>((**std::shared_ptr<std::shared_ptr<Symbol> >(buf_[i])).getRef());
            // size estimate, not guaranteed to be totally accurate
            // because of the weak_ptr contribution
            (*mem) += static_cast<size_t>(
                    (static_cast<long double>(len * sizeof(char) + sizeof(Symbol))
                     * one_over_secondary_ref_count)
                    + static_cast<long double>(sizeof(std::weak_ptr<std::shared_ptr<Symbol> >)) );
        }
    }
}

void Symbuf::getTopology(std::map<long int, std::map<std::string, std::map<std::shared_ptr<Symbol>*, size_t> > > &topology) const
{

    // topology: group ID ->
    //                string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                          -> reference counts of the shared pointer sets
    typedef std::shared_ptr<Symbol>*         optr_type;
    typedef std::map<optr_type, size_t>      imap_type;
    typedef std::map<std::string, imap_type> mmap_type;
    typedef std::map<long int, mmap_type>    omap_type;

    for (size_t i = 0; i < bufsize_; ++i)
    {
        if (buf_[i].expired())
        {
            continue;
        }
        else
        {
            std::shared_ptr<std::shared_ptr<Symbol> > tmp_ptr(buf_[i]);
            std::string         this_string((**tmp_ptr).getPtr());
            long int            this_group((**tmp_ptr).getGroup());

            omap_type::iterator it1 = topology.find(this_group);
            // look for an existing entry for this group
            if (it1 == topology.end())
            {
                // create a new group entry
                mmap_type group_tmp_map;
                it1 = topology.insert(topology.end(), omap_type::value_type(this_group, group_tmp_map));
            }

            // look for an existing entry for this_string
            mmap_type::iterator it2 = it1->second.find(this_string);
            if (it2 == it1->second.end())
            {
                // entry not yet in map -> create a new entry
                imap_type tmp_map;
                tmp_map.insert(imap_type::value_type(tmp_ptr.get(), (tmp_ptr.use_count() - 1)));
                it1->second.insert(mmap_type::value_type(this_string, tmp_map));
            }
            else
            {
                // the entry already exists -> expand it
                // look for an existing entry for tmp_ptr
                imap_type::iterator it3 = it2->second.find(tmp_ptr.get());
                if (it3 != it2->second.end())
                {
                    // shared ptr set is already present, but the reference count may have changed in the meantime
                    it3->second = tmp_ptr.use_count();
                }
                else
                {
                    // this shared_ptr set is not yet registered, create a new entry
                    it2->second.insert(imap_type::value_type(tmp_ptr.get(), (tmp_ptr.use_count() - 1)));
                }
            }
        }
    }

}

// Symtab member functions

// trim white space at the string ends and remaining chars at the end exceeding maxlen
std::string Symtab::trimString(const char * const s1, const size_t maxlen) const
{
    if (s1 == nullptr)
    {
        return std::string("");
    }

    int  i = 0,  ii = 0;

    // skip leading whitespace
    while (&s1[i] && (s1[i] == ' ' || s1[i] == '\t'))
    {
        i++;
    }
    // find string end (position of the last non-terminating characters)
    ii = i;
    while (&s1[ii+1] != nullptr && s1[ii+1] != '\0')
    {
        ii++;
    }
    // omit trailing whitespace
    while (ii > i && (s1[ii] == ' ' || s1[ii] == '\t'))
    {
        ii--;
    }

    size_t len = ii - i + 1;
    if (len >= maxlen)
    {
        len = maxlen - 1;
        fprintf(stderr, "Warning in %s::trim_string():\n", typeid(*this).name());
        fprintf(stderr, "Truncated exceedingly long std::string to max. length %li, std::string: %s\n", maxlen, s1);
    }
    // copy the trimmed string to a char array
    char * const tmp_str = new char[len+1];
    const char * s       = &s1[i];
    for (i = 0; i < (int)len; i++)
    {
        tmp_str[i] = *(s++);
    }
    tmp_str[i] = '\0';

    std::string result(tmp_str);
    delete [] tmp_str;

    return result;
}

/*
   utility function for symbol table, trims whitespace at the beginning of
   a string and trims the string to a total maximum length given by maxlen
   the resulting string is written to an output string that is managed by a shared_ptr
   with the correct size
 */
std::shared_ptr<std::shared_ptr<Symbol> > Symtab::constructEntry(const char * const in, const long int group_in, const size_t maxlen, size_t *tmp_mem)
{
    if (in == nullptr)
    {
        return nullptr;
    }

    size_t      len, i;
    std::string tmp_str(trimString(in, maxlen));
    const char *s = tmp_str.c_str();

    while ((*s) && ((*s) == ' ' || (*s) == '\t'))
    {
        s++;
    }
    for (len = std::strlen(s); (len > 0); len--)
    {
        if (s[len] != ' ')
        {
            break;
        }
    }
    if (len >= maxlen)
    {
        len = maxlen - 1;
        fprintf(stderr, "Warning in  %s::constructEntry()\n", typeid(*this).name());
        fprintf(stderr, "Truncated exceedingly long std::string to max. length %lu. std::string: %s\n", static_cast<unsigned long int>(maxlen), in);
        fprintf(stderr, "Error: Insufficient memory for allocating new symbol table buffer!\n");
    }
    // functions for memory allocation/deallocation of the shared_ptr
    // ... for the object pointed to and managed by the outer pointer, i.e.,
    //     for creating and disposing of the inner pointer
    // ... for the object pointed to and managed by the inner pointer
    auto innerAllocator = [&]() -> Symbol* {
            Symbol* p = new Symbol;
            p->setRef(1);
            /* p->setPtr(new char[len+1]); */
            char *cp = (new char[len+1]);
            cp[0] = '\0';
            p->setPtr(cp);
            p->setGroup(group_in);
            return p;
        };
    // safe variant (?)
    // 1. atomically swap the address of the allocated memory with a temporary pointer that holds nullptr
    //    -> allocated memory is no longer accessible to other threads
    // 2. deallocate the memory via delete [] acting on the temporary pointer
    // 3. delete the Symbol enclosing the atomic pointer itself
    auto innerDeleter   = [&](Symbol* p) -> void {
            if (p != nullptr)
            {
                if (p->getRef() == 0)
                {
                    const char* tmp_p = p->exchangePtr(nullptr);
                    if (tmp_p != nullptr)
                    {
                        delete [] tmp_p;
                    }
                    delete p;
                }
            }
        };
    auto outerAllocator = [&]() -> std::shared_ptr<Symbol>* {
            std::shared_ptr<Symbol> *p =
                new std::shared_ptr<Symbol>( innerAllocator(), innerDeleter );
            return p;
        };
    auto outerDeleter   = [&](std::shared_ptr<Symbol>* p) -> void {
            if (p != nullptr)
            {
                // decrement the secondary reference counter,
                // if the resulting secondary reference count
                // is zero, no other outer pointer referencing
                // the object is present and we can delete it
                if ((*p) != nullptr)
                {
                    size_t rcnt = (*p)->decrRef();
                    if (rcnt == 0)
                    {
                    }
                }
                delete p;
            }
        };
    std::shared_ptr<std::shared_ptr<Symbol> > result( outerAllocator(), outerDeleter );

    char * tmp_cp = (*result)->getPtr();
    // store the chars forming the string
    for (i = 0; i < len; i++)
    {
        tmp_cp[i] = *(s++);
    }
    tmp_cp[i] = '\0';
    // memory occupancy
    (*tmp_mem) = (len + 1) * sizeof(char) + sizeof(Symbol);
    return result;
}

// the function does not lock the this->mtx_ by itself because it is called
// by other functions that should already own the lock
std::shared_ptr<std::shared_ptr<Symbol> > Symtab::nonconstGetEntry(const char *s, const long int group) const
{
    if (s == nullptr)
    {
        return nullptr;
    }

    // all entries with group IDs < 0 belong to the same default string group
    long int gid = group;
    if (gid < 0)
    {
        gid = -1;
    }

    std::shared_ptr<Symbuf>                   search_pt = symbuf_;
    std::shared_ptr<std::shared_ptr<Symbol> > result;

    while (search_pt != nullptr)
    {
        result = search_pt->getEntry(s, gid);
        if (result != nullptr)
        {
            return result;
        }
        search_pt = search_pt->next_;
    }

    return result;
}

std::shared_ptr<std::shared_ptr<Symbol> > Symtab::nonconstInsertEntry(const char * const s, const long int group)
{
    if (s == nullptr)
    {
        return nullptr;
    }

    // make sure that all entries with group IDs < 0 are added to the same default string group
    long int gid = group;
    if (gid < 0)
    {
        gid = -1;
    }

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    std::lock_guard<std::mutex> lock(mtx_);

    // first try if the entry exists already, otherwise create a new entry
    std::shared_ptr<std::shared_ptr<Symbol> > result(nonconstGetEntry(s, gid));
    if (result != nullptr)
    {
        return result;
    }
    else
    {

        // look for an existing buffer with an unoccupied entry
        // and if it exists store the new entry there
        std::shared_ptr<Symbuf> this_pt = symbuf_;
        std::shared_ptr<Symbuf> prev_pt = nullptr;
        size_t                  tmp_mem;

        result = constructEntry(s, gid, maxEntrySize_, &tmp_mem);
        bool success = false;
        while (this_pt != nullptr)
        {
            try
            {
                success = this_pt->insertEntry(result);
            }
            catch (std::logic_error &)
            {
                fprintf(stderr, "[exception caught] in %s::insertEntry():\n", typeid(*this).name());
                fprintf(stderr, "failed to insert entry because of internal inconsistency\n");
                throw;
            }
            catch (std::bad_alloc &)
            {
                fprintf(stderr, "[exception caught] in %s::insertEntry():\n", typeid(*this).name());
                fprintf(stderr, "failed to allocate memory for std::string entry\n");
                throw;
            }

            if (success)
            {
                nr_++;
                mem_ += tmp_mem;
                return result;
            }
            prev_pt = this_pt;
            this_pt = this_pt->next_;
        }

        // allocate memory for a new buffer because
        // all preexisting buffers were full
        if (this_pt == nullptr)
        {
            this_pt = prev_pt;
            // we could still have no Symbuf at all
            if (this_pt == nullptr)
            {
                symbuf_  = std::make_shared<Symbuf>();
                this_pt  = symbuf_;
            }
            else
            {
                this_pt->next_ = std::make_shared<Symbuf>();
                this_pt        = this_pt->next_;
            }
            this_pt->next_ = nullptr;
            this_pt->prev_ = prev_pt;

            try
            {
                this_pt->bufsize_ = bufsize_;
                tmp_mem           = this_pt->allocSymbuf();
            }
            catch (std::bad_alloc &)
            {
                fprintf(stderr, "Error in  %s::insertEntry()\n", typeid(*this).name());
                fprintf(stderr, "Insufficient memory for allocating new symbol table buffer!\n");
                throw;
            }

            nbuf_++;
            mem_ += tmp_mem;

        }
        else
        {
            fprintf(stderr, "Error in  %s::insertEntry()\n", typeid(*this).name());
            fprintf(stderr, "Failed to find end of linked list for allocating a new symbol table buffer for unknown reason.\n");
            throw std::logic_error("Internal error in Symtab");
        }

        // store the new string entry
        success = this_pt->insertEntry(result);

        nr_++;
        mem_ += tmp_mem;

        // check whether the entry was inserted successfully
        if (!success)
        {
            // this should never happen as long as the code of Symtab and
            // Symbuf is internally consistent
            fprintf(stderr, "Error in  %s::insertEntry()\n", typeid(*this).name());
            fprintf(stderr, "Could not insert new entry for unknown reason.\n");
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
        return result;
    }

//// if (result) {
////   std::cout << "  1st level is OK. use count " << result.useCount() << std::endl;
////   if (*result) {
////     std::cout << "    2nd level is OK. use count " << result.useCount() << std::endl;
////     std::cout << "    (*result)->getRef()       : " << (*result)->getRef() << std::endl;
////     std::cout << "    (*result)->getPtr()       : " << static_cast<void*>((*result)->getPtr()) << std::endl;
////     if ((*result)->getPtr()) {
////       std::cout << "    (*result)->getPtr() string: " << (*result)->getPtr() << std::endl;
////     }
////   } else {
////     std::cout << "    2nd level is NOT OK." << std::endl;
////   }
//// } else {
////   std::cout << "  1st level is NOT OK." << std::endl;
//// }

    return nullptr;
}

std::shared_ptr<std::shared_ptr<Symbol> > Symtab::nonconstRenameEntry(const char * const before, const long int group_before,
                                                                      const char * const after,  const long int group_after)
{
    if (before == nullptr || after == nullptr)
    {
        return nullptr;
    }

    // make sure that all entries with group IDs < 0 are added to the same default group
    long int gid_before = group_before;
    if (gid_before < 0)
    {
        gid_before = -1;
    }
    long int gid_after = group_after;
    if (gid_after < 0)
    {
        gid_after = -1;
    }

    // the nesting of multiple member functions necessitates manual locking
    // and unlocking instead of the convenient use of lock_guard
    // changes by other threads could take place in between the periods of
    // owning the lock but that should not cause problems as long as there
    // are no simultaneous calls to the function by multiple threads
    // first try whether the entry before exists already,
    // otherwise delegate the job to insert_entry()
    mtx_.lock();
    std::shared_ptr<std::shared_ptr<Symbol> > result(nonconstGetEntry(before, gid_before));
    mtx_.unlock();

    if (result != nullptr)
    {

        // If the entry "after" already exists, exchange the target of the
        // pointer pointing to "before". Otherwise, store the new string
        // "after" in a newly allocated memory area and let the pointer
        // to "before" point to it.
        std::shared_ptr<std::shared_ptr<Symbol> > tmp_ptr(nonconstInsertEntry(after, gid_after));
        if (tmp_ptr != nullptr)
        {

            // other threads might have made changes that removed the
            // previously found entry for before in the meantime
            // "before" is preserved by the referencing shared_ptr instance "result"
            // "after", by "tmp_ptr"
            std::lock_guard<std::mutex> lock(this->mtx_);
            result = nonconstGetEntry(before, gid_before);
            while (result != nullptr)
            {
                // There seems to be no way to manipulate the control block of the shared pointers
                // such that two sets of shared_pointer instances can be merged, a custom shared
                // pointer with that capability should in principle be possible but detach the code
                // from possible future improvements of std::shared_ptr. Is there a safe, strictly
                // portable way via inheritance from shared_ptr? probably not.
                (*tmp_ptr)->incrRef();
                std::shared_ptr<Symbol> save_last(*result);
                *result = *tmp_ptr;
                result  = nonconstGetEntry(before, gid_before);
                if (result == nullptr)
                {
                    save_last->setRef((size_t)0);
                }
            }

            result = nonconstGetEntry(after, gid_after);
            return result;
        }
        else
        {
            // this should never happen as long as the code of Symtab and
            // Symbuf is internally consistent
            fprintf(stderr, "Error in  %s::nonconstRenameEntry()\n", typeid(*this).name());
            fprintf(stderr, "Could not insert or locate entry for unknown reason.\n");
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
    }
    else
    {
        // not a fatal error, simply create a new entry or return a shared_ptr
        // to a possible preexisting entry for "after"
        return nonconstInsertEntry(after, gid_after);
    }

    return nullptr;
}

void Symtab::getStringStats(const char * const s, const long int group, size_t *ninst, size_t *nptr) const
{
    if (s == nullptr)
    {
        return;
    }

    // make sure that all entries with group IDs < 0 refer to the same default group
    long int gid = group;
    if (gid < 0)
    {
        gid = -1;
    }

    // topology: group ID ->
    //                string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                          -> reference counts of the shared pointer sets
    typedef std::shared_ptr<Symbol>*         optr_type;
    typedef std::map<optr_type, size_t>      imap_type;
    typedef std::map<std::string, imap_type> mmap_type;
    typedef std::map<long int, mmap_type>    omap_type;
    omap_type           topology;
    getTopology(topology);

    omap_type::iterator it1 = topology.find(gid);
    if (it1 != topology.end())
    {
        mmap_type::iterator it2 = it1->second.find(std::string(s));

        *nptr  = 0;
        *ninst = 0;
        if (it2 != it1->second.end())
        {

            for (imap_type::iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
            {
                (*nptr)++;
                (*ninst) += it3->second;
            }

            if (*nptr != (*nonconstGetEntry(s, gid))->getRef())
            {
                // this should never happen as long as the code of Symtab and
                // Symbuf is internally consistent
                fprintf(stderr, "Error in  %s::getStringStats()\n", typeid(*this).name());
                fprintf(stderr, "Detected inconsistent secondary reference counts.\n");
                fprintf(stderr, " getTopology(): %li vs. Symbol: %li\n", *nptr, (*nonconstGetEntry(s))->getRef());
            }

        }
        else
        {
            // this should never happen as long as the code of Symtab and
            // Symbuf is internally consistent
            fprintf(stderr, "Error in  %s::getStringStats()\n", typeid(*this).name());
            fprintf(stderr, "Could not find entry: \"%s\"\n", s);
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
    }
    else
    {
        // this should never happen as long as the code of Symtab and
        // Symbuf is internally consistent
        fprintf(stderr, "Error in  %s::getStringStats()\n", typeid(*this).name());
        fprintf(stderr, "Could not find entry group: \"%s\"\n", s);
        throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
    }

}

void Symtab::getTopology(std::map<long int, std::map<std::string, std::map<std::shared_ptr<Symbol>*, size_t> > > &topology) const
{

    std::shared_ptr<Symbuf> search_pt = symbuf_;

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    std::lock_guard<std::mutex> lock(this->mtx_);

    while (search_pt != nullptr)
    {
        search_pt->getTopology(topology);
        search_pt = search_pt->next_;
    }

}

void Symtab::updateStats()
{

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    std::lock_guard<std::mutex> lock(mtx_);
    mem_  = sizeof(Symtab);
    nr_   = 0;
    nbuf_ = 0;

    // only proceed if there is any symbol buffer yet
    if (symbuf_ == nullptr)
    {
        return;
    }

    // sum up the numbers of Symbufs and their string entries
    // and the memory occupied by them
    size_t                    tmp_mem = 0;
    unsigned int              tmp_nr  = 0;
    std::shared_ptr<Symbuf>   tmp_buf1(symbuf_);
    while (tmp_buf1 != nullptr)
    {
        tmp_buf1->getStats(&tmp_nr, &tmp_mem);
        nr_  += tmp_nr;
        mem_ += tmp_mem;
        nbuf_++;
        tmp_buf1 = tmp_buf1->next_;
    }
}

void Symtab::clean()
{

    // return if there is no stored data at all
    if (symbuf_ == nullptr)
    {
        printf("Symtab::cleanSymtab(): returning because there is no data allocated\n");
        return;
    }

    // determine number of entries and symbol buffers
    updateStats();

    //size_t tmp_mem = mem_;

    this->mtx_.lock();

    std::shared_ptr<Symbuf>   tmp_buf1(symbuf_);
    size_t                    nmoved   = 0;
    size_t                    nr_avail = 0;
    while (tmp_buf1->next_ != nullptr)
    {
        tmp_buf1 = tmp_buf1->next_;
        nbuf_++;
        nr_avail += tmp_buf1->bufsize_;
    }
    // return if all buffers are completely filled
    if (nr_avail == nr_)
    {
        printf("Symtab::cleanSymtab(): returning because all buffers are completely filled and hence unfragmented\n");
        return;
    }

    bool fragmented = true;
    while (fragmented)
    {
        fragmented = false;

        // find the last occupied position
        size_t                    tmp_index1   = 0;
        bool                      got_occupied = false;
        std::shared_ptr<Symbuf>   last_occupied(symbuf_);
        std::shared_ptr<Symbuf>   ibuf(symbuf_);
        while (ibuf->next_ != nullptr)
        {
            ibuf = ibuf->next_;
            for (size_t i = 0; i < ibuf->bufsize_; i++)
            {
                if (!ibuf->buf_[i].expired())
                {
                    last_occupied = ibuf;
                    tmp_index1    = i;
                    got_occupied  = true;
                }
            }
        }

        if (!got_occupied)
        {
            break;
        }

        // find the first vacant position
        size_t                  tmp_index2 = 0;
        std::shared_ptr<Symbuf> first_free(symbuf_);
        ibuf = symbuf_;
        while (ibuf->next_ != nullptr)
        {
            ibuf = ibuf->next_;
            for (size_t i = 0; i < ibuf->bufsize_; i++)
            {
                // stop at the previously determined last occupied position
                if (i == tmp_index1 && ibuf == last_occupied)
                {
                    break;
                }

                if (ibuf->buf_[i].expired())
                {
                    first_free = ibuf;
                    tmp_index2 = i;
                    fragmented = true;
                    break;
                }
            }
        }

        if (fragmented)
        {
            first_free->buf_[tmp_index2].swap(last_occupied->buf_[tmp_index1]);
            nmoved++;
        }

    }   // ... end while (fragmented) {}

    // delete empty buffers beginning at the end of the liked list
    tmp_buf1 = symbuf_;
    size_t nbufdeleted = 0;
    while (tmp_buf1->next_ != nullptr)
    {
        tmp_buf1 = tmp_buf1->next_;
    }
    while (tmp_buf1 != symbuf_)
    {
        if (tmp_buf1->buf_ == nullptr)
        {
            break;
        }
        else if (tmp_buf1->buf_[0].expired())
        {
            break;
        }
        tmp_buf1 = tmp_buf1->prev_.lock();
        tmp_buf1->next_.reset();
        nbufdeleted++;
    }

    this->mtx_.unlock();

    // update current memory occupied, number of symbol buffers
    // and number of std::string entries
    updateStats();

/*
    // instead of commenting out to avoid compiler warnings because of unused variables
    tmp_mem = mem_ - tmp_mem;
    if (1 == 0)
    {
        printf("%s::clean():\n", typeid(*this).name());
        printf("moved %li entries, and deleted %li buffers\n", nmoved, nbufdeleted);
        printf("freed %li bytes of memory.\n", tmp_mem);
    }
 */
}

} // namespace gmx
