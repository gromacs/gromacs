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
    Implementation file needed for linking all code parts that share
    symbol table entries to a common symbol table object.

   \author
    R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \copyright
    GROMACS license

   \date Sep 2015
 */

#include "gromacs/data_structures/symtab.h"

namespace gmx
{

#if    defined gmx_symtab_cpp11_h &&  defined gmx_symtab_cpp98_h
 #error "Headers for both symbol table variants were included"
#elif !defined gmx_symtab_cpp11_h && !defined gmx_symtab_cpp98_h
 #error "No header for any symbol table variant was included"
#elif  defined gmx_symtab_cpp11_h

//=====================================================================================
//
// implementations of member functions for the C++11 variant
//
//=====================================================================================

/*! function compares two strings, not considering wildcards, and
    returns true if the strings are equal and false otherwise
    leading and trailing whitespace is skipped                    */
bool t_symbuf::match_strings(const char * const s1, const char * const s2) const
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

bool t_symbuf::insert_entry(std::shared_ptr<std::shared_ptr<t_symbol> > in)
{
    // Search for a free pointer in the array buf.
    // If found, allocate memory for the new string entry
    // and store it in the newly allocated space
    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i].expired())
        {
            buf[i] = in;
            return true;
        }
    }

    return false;
}

std::shared_ptr<std::shared_ptr<t_symbol> > t_symbuf::get_entry(const char * const s) const
{
    if (s == nullptr)
    {
        return nullptr;
    }

    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i].expired())
        {
            continue;
        }
        else
        {
            std::shared_ptr<std::shared_ptr<t_symbol> > tmp_ptr(buf[i]);
            if (match_strings((**tmp_ptr).get_ptr(), s))
            {
                return tmp_ptr;
            }
        }
    }

    return nullptr;
}

void t_symbuf::get_stats(unsigned int * const nr, size_t * const mem)
{
    (*nr)  = 0;
    (*mem) = sizeof(t_symbuf) + sizeof(std::weak_ptr<std::shared_ptr<t_symbol> >) * bufsize;
    for (size_t i = 0; i < bufsize; ++i)
    {
        if (!buf[i].expired())
        {
            (*nr)++;
            size_t      len = string_length((**std::shared_ptr<std::shared_ptr<t_symbol> >(buf[i])).get_ptr());
            long double one_over_secondary_ref_count = 1 /
                static_cast<long double>((**std::shared_ptr<std::shared_ptr<t_symbol> >(buf[i])).get_ref());
            // size estimate, not guaranteed to be totally accurate
            // because of the weak_ptr contribution
            (*mem) += static_cast<size_t>(
                    (static_cast<long double>(len * sizeof(char) + sizeof(t_symbol))
                     * one_over_secondary_ref_count)
                    + static_cast<long double>(sizeof(std::weak_ptr<std::shared_ptr<t_symbol> >)) );
        }
    }
}

void t_symbuf::get_topology(std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> > &topology) const
{

    // topology: string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
    typedef std::shared_ptr<t_symbol>* optr;

    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i].expired())
        {
            continue;
        }
        else
        {
            std::shared_ptr<std::shared_ptr<t_symbol> > tmp_ptr(buf[i]);
            std::string this_string((**tmp_ptr).get_ptr());
            // look for an existing entry for this_string
            std::map<std::string, std::map<optr, size_t> >::iterator it1 = topology.find(this_string);
            if (it1 == topology.end())
            {
                // entry not yet in map -> create a new entry
                std::map<optr, size_t> tmp_map;
                tmp_map.insert(std::map<optr, size_t>::value_type(tmp_ptr.get(), (tmp_ptr.use_count() - 1)));
                topology.insert(std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> >::value_type(this_string, tmp_map));
            }
            else
            {
                // the entry already exists -> expand it
                // look for an existing entry for tmp_ptr
                std::map<optr, size_t>::iterator it2 = it1->second.find(tmp_ptr.get());
                if (it2 != it1->second.end())
                {
                    // shared ptr set is already present, but the reference count may have changed in the meantime
                    it2->second = tmp_ptr.use_count();
                }
                else
                {
                    // this shared_ptr set is not yet registered, create a new entry
                    it1->second.insert(std::map<optr, size_t>::value_type(tmp_ptr.get(), (tmp_ptr.use_count() - 1)));
                }
            }
        }
    }

}

// t_symtab member functions

//! trim white space at the string ends and remaining chars at the end exceeding maxlen
std::string t_symtab::trim_string(const char * const s1, const size_t maxlen) const
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

std::shared_ptr<std::shared_ptr<t_symbol> > t_symtab::construct_entry(const char * const in, const size_t maxlen, size_t *tmp_mem)
/*!
 * utility function for symbol table, trims whitespace at the beginning of
 * a string and trims the string to a total maximum length given by maxlen
 * the resulting string is written to an output string that is managed by a shared_ptr
 * with the correct size
 */
{
    if (in == nullptr)
    {
        return nullptr;
    }

    size_t      len, i;
    std::string tmp_str(trim_string(in, maxlen));
    const char *s = tmp_str.c_str();

    while ((*s) && ((*s) == ' ' || (*s) == '\t'))
    {
        s++;
    }
    for (len = string_length(s) - 1; (len > 0); len--)
    {
        if (s[len] != ' ')
        {
            break;
        }
    }
    if (len >= maxlen)
    {
        len = maxlen - 1;
        fprintf(stderr, "Warning in  %s::construct_entry()\n", typeid(*this).name());
        fprintf(stderr, "Truncated exceedingly long std::string to max. length %li. std::string: %s\n", maxlen, in);
        fprintf(stderr, "Error: Insufficient memory for allocating new symbol table buffer!\n");
    }
    // functions for memory allocation/deallocation of the shared_ptr
    // ... for the object pointed to and managed by the outer pointer, i.e.,
    //     for creating and disposing of the inner pointer
    // ... for the object pointed to and managed by the inner pointer
    auto inner_allocator = [&]() -> t_symbol* {
            t_symbol* p = new t_symbol;
            p->set_ref(1);
            /* p->set_ptr(new char[len+1]); */
            char *cp = (new char[len+1]);
            cp[0] = '\0';
            p->set_ptr(cp);
            return p;
        };
    // safe variant (?)
    // 1. atomically swap the address of the allocated memory with a temporary pointer that holds nullptr/NULL
    //    -> allocated memory is no longer accessible to other threads
    // 2. deallocate the memory via delete [] acting on the temporary pointer
    // 3. delete the t_symbol enclosing the atomic pointer itself
    auto inner_deleter   = [&](t_symbol* p) -> void {
            if (p != nullptr)
            {
                if (p->get_ptr() != nullptr)
                {
                }
                else
                {
                }
                if (p->get_ref() == 0)
                {
                    const char* tmp_p = p->exchange_ptr(nullptr);
                    if (tmp_p != nullptr)
                    {
                        delete [] tmp_p;
                    }
                    delete p;
                }
            }
        };
    auto outer_allocator = [&]() -> std::shared_ptr<t_symbol>* {
            std::shared_ptr<t_symbol> *p =
                new std::shared_ptr<t_symbol>( inner_allocator(), inner_deleter );
            return p;
        };
    auto outer_deleter   = [&](std::shared_ptr<t_symbol>* p) -> void {
            if (p != nullptr)
            {
                // decrement the secondary reference counter,
                // if the resulting secondary reference count
                // is zero, no other outer pointer referencing
                // the object is present and we can delete it
                if ((*p) != nullptr)
                {
                    size_t rcnt = (*p)->decr_ref();
                    if (rcnt == 0)
                    {
                    }
                }
                delete p;
            }
        };
    std::shared_ptr<std::shared_ptr<t_symbol> > result( outer_allocator(), outer_deleter );

    char * tmp_cp = (*result)->get_ptr();
    // store the chars forming the string
    for (i = 0; i < len; i++)
    {
        tmp_cp[i] = *(s++);
    }
    tmp_cp[i] = '\0';
    // memory occupancy
    (*tmp_mem) = (len + 1) * sizeof(char) + sizeof(t_symbol);
    return result;
}

// the function does not lock the this->mtx by itself because it is called
// by other functions that should already own the lock
std::shared_ptr<std::shared_ptr<t_symbol> > t_symtab::nonconst_get_entry(const char *s) const
{
    if (s == nullptr)
    {
        return nullptr;
    }

    std::shared_ptr<t_symbuf>                   search_pt = symbuf;
    std::shared_ptr<std::shared_ptr<t_symbol> > result;

    while (search_pt != nullptr)
    {
        result = search_pt->get_entry(s);
        if (result != nullptr)
        {
            return result;
        }
        search_pt = search_pt->next;
    }

    return result;
}

std::shared_ptr<std::shared_ptr<t_symbol> > t_symtab::nonconst_insert_entry(const char * const s)
{
    if (s == nullptr)
    {
        return nullptr;
    }

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    std::lock_guard<std::mutex> lock(this->mtx);

    // first try if the entry exists already, otherwise create a new entry
    std::shared_ptr<std::shared_ptr<t_symbol> > result(nonconst_get_entry(s));
    if (result != nullptr)
    {
        return result;
    }
    else
    {

        // look for an existing buffer with an unoccupied entry
        // and if it exists store the new entry there
        std::shared_ptr<t_symbuf> this_pt = symbuf;
        std::shared_ptr<t_symbuf> prev_pt = nullptr;
        size_t                    tmp_mem;

        result = construct_entry(s, max_entry_size, &tmp_mem);
        bool success = false;
        while (this_pt != nullptr)
        {
            try
            {
                success = this_pt->insert_entry(result);
            }
            catch (std::logic_error &)
            {
                fprintf(stderr, "[exception caught] in %s::insert_entry():\n", typeid(*this).name());
                fprintf(stderr, "failed to insert entry because of internal inconsistency\n");
                throw;
            }
            catch (std::bad_alloc &)
            {
                fprintf(stderr, "[exception caught] in %s::insert_entry():\n", typeid(*this).name());
                fprintf(stderr, "failed to allocate memory for std::string entry\n");
                throw;
            }

            if (success)
            {
                nr++;
                mem += tmp_mem;
                return result;
            }
            prev_pt = this_pt;
            this_pt = this_pt->next;
        }

        // allocate memory for a new buffer because
        // all preexisting buffers were full
        if (this_pt == nullptr)
        {
            this_pt = prev_pt;
            // we could still have no t_symbuf at all
            if (this_pt == nullptr)
            {
                symbuf  = std::make_shared<t_symbuf>();
                this_pt = symbuf;
            }
            else
            {
                this_pt->next = std::make_shared<t_symbuf>();
                this_pt       = this_pt->next;
            }
            this_pt->next = nullptr;
            this_pt->prev = prev_pt;

            try
            {
                this_pt->bufsize = bufsize;
                this_pt->alloc_symbuf(&tmp_mem);
            }
            catch (std::bad_alloc &)
            {
                fprintf(stderr, "Error in  %s::insert_entry()\n", typeid(*this).name());
                fprintf(stderr, "Insufficient memory for allocating new symbol table buffer!\n");
                throw;
            }
            ;
            nbuf++;
            mem += tmp_mem;

        }
        else
        {
            fprintf(stderr, "Error in  %s::insert_entry()\n", typeid(*this).name());
            fprintf(stderr, "Failed to find end of linked list for allocating a new symbol table buffer for unknown reason.\n");
            throw std::logic_error("Internal error in t_symtab");
        }

        // store the new string entry
        success = this_pt->insert_entry(result);

        nr++;
        mem += tmp_mem;

        // check whether the entry was inserted successfully
        if (!success)
        {
            // this should never happen as long as the code of t_symtab and
            // t_symbuf is internally consistent
            fprintf(stderr, "Error in  %s::insert_entry()\n", typeid(*this).name());
            fprintf(stderr, "Could not insert new entry for unknown reason.\n");
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
        return result;
    }

//// if (result) {
////   std::cout << "  1st level is OK. use count " << result.use_count() << std::endl;
////   if (*result) {
////     std::cout << "    2nd level is OK. use count " << result.use_count() << std::endl;
////     std::cout << "    (*result)->get_ref()       : " << (*result)->get_ref() << std::endl;
////     std::cout << "    (*result)->get_ptr()       : " << static_cast<void*>((*result)->get_ptr()) << std::endl;
////     if ((*result)->get_ptr()) {
////       std::cout << "    (*result)->get_ptr() string: " << (*result)->get_ptr() << std::endl;
////     }
////   } else {
////     std::cout << "    2nd level is NOT OK." << std::endl;
////   }
//// } else {
////   std::cout << "  1st level is NOT OK." << std::endl;
//// }

    return nullptr;
}

std::shared_ptr<std::shared_ptr<t_symbol> > t_symtab::nonconst_rename_entry(const char * const before, const char * const after)
{
    if (before == nullptr || after == nullptr)
    {
        return nullptr;
    }

    // the nesting of multiple member functions necessitates manual locking
    // and unlocking instead of the convenient use of lock_guard
    // changes by other threads could take place in between the periods of
    // owning the lock but that should not cause problems as long as there
    // are no simultaneous calls to the function by multiple threads
    // first try whether the entry before exists already,
    // otherwise delegate the job to insert_entry()
    this->mtx.lock();
    std::shared_ptr<std::shared_ptr<t_symbol> > result(nonconst_get_entry(before));
    this->mtx.unlock();

    if (result != nullptr)
    {

        // If the entry "after" already exists, exchange the target of the
        // pointer pointing to "before". Otherwise, store the new string
        // "after" in a newly allocated memory area and let the pointer
        // to "before" point to it.
        std::shared_ptr<std::shared_ptr<t_symbol> > tmp_ptr(nonconst_insert_entry(after));
        if (tmp_ptr != nullptr)
        {

            // other threads might have made changes that removed the
            // previously found entry for before in the meantime
            // "before" is preserved by the referencing shared_ptr instance "result"
            // "after", by "tmp_ptr"
            std::lock_guard<std::mutex> lock(this->mtx);
            result = nonconst_get_entry(before);
            while (result != nullptr)
            {
                // There seems to be no way to manipulate the control block of the shared pointers
                // such that two sets of shared_pointer instances can be merged, a custom shared
                // pointer with that capability should in principle be possible but detach the code
                // from possible future improvements of std::shared_ptr. Is there a safe, strictly
                // portable way via inheritance from shared_ptr? probably not.
                (*tmp_ptr)->incr_ref();
                std::shared_ptr<t_symbol> save_last(*result);
                *result = *tmp_ptr;
                result  = nonconst_get_entry(before);
                if (result == nullptr)
                {
                    save_last->set_ref((size_t)0);
                }
            }

            result = nonconst_get_entry(after);
            return result;
        }
        else
        {
            // this should never happen as long as the code of t_symtab and
            // t_symbuf is internally consistent
            fprintf(stderr, "Error in  %s::nonconst_rename_entry()\n", typeid(*this).name());
            fprintf(stderr, "Could not insert or locate entry for unknown reason.\n");
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
    }
    else
    {
        // not a fatal error, simply create a new entry or return a shared_ptr
        // to a possible preexisting entry for "after"
        return nonconst_insert_entry(after);
    }

    return nullptr;
}

void t_symtab::get_string_stats(const char * const s, size_t *ninst, size_t *nptr) const
{
    if (s == nullptr)
    {
        return;
    }

    // topology: string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
    typedef std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> > map_type;
    typedef std::map<std::shared_ptr<t_symbol>*, size_t> imap_type;
    map_type           topology;
    get_topology(topology);
    map_type::iterator it1 = topology.find(std::string(s));

    *nptr  = 0;
    *ninst = 0;
    if (it1 != topology.end())
    {

        for (imap_type::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
            (*nptr)++;
            (*ninst) += it2->second;
        }

        if (*nptr != (*nonconst_get_entry(s))->get_ref())
        {
            // this should never happen as long as the code of t_symtab and
            // t_symbuf is internally consistent
            fprintf(stderr, "Error in  %s::get_string_stats()\n", typeid(*this).name());
            fprintf(stderr, "Detected inconsistent secondary reference counts.\n");
            fprintf(stderr, " get_topology(): %li vs. t_symbol: %li\n", *nptr, (*nonconst_get_entry(s))->get_ref());
        }

    }
    else
    {
        // this should never happen as long as the code of t_symtab and
        // t_symbuf is internally consistent
        fprintf(stderr, "Error in  %s::get_string_stats()\n", typeid(*this).name());
        fprintf(stderr, "Could not find entry: \"%s\"\n", s);
        throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
    }

}

void t_symtab::get_topology(std::map<std::string, std::map<std::shared_ptr<t_symbol>*, size_t> > &topology) const
{

    std::shared_ptr<t_symbuf> search_pt = symbuf;

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    std::lock_guard<std::mutex> lock(this->mtx);

    while (search_pt != nullptr)
    {
        search_pt->get_topology(topology);
        search_pt = search_pt->next;
    }

}

void t_symtab::update_stats()
{

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    std::lock_guard<std::mutex> lock(this->mtx);
    mem  = sizeof(t_symtab);
    nr   = 0;
    nbuf = 0;

    // only proceed if there is any symbol buffer yet
    if (symbuf == nullptr)
    {
        return;
    }

    // sum up the numbers of t_symbufs and their string entries
    // and the memory occupied by them
    size_t                    tmp_mem = 0;
    unsigned int              tmp_nr  = 0;
    std::shared_ptr<t_symbuf> tmp_buf1(symbuf);
    while (tmp_buf1 != nullptr)
    {
        tmp_buf1->get_stats(&tmp_nr, &tmp_mem);
        nr  += tmp_nr;
        mem += tmp_mem;
        nbuf++;
        tmp_buf1 = tmp_buf1->next;
    }
}

void t_symtab::clean()
{

    // return if there is no stored data at all
    if (symbuf == NULL)
    {
        printf("t_symtab::clean_symtab(): returning because there is no data allocated\n");
        return;
    }

    // determine number of entries and symbol buffers
    update_stats();
    size_t tmp_mem = mem;

    this->mtx.lock();

    std::shared_ptr<t_symbuf>   tmp_buf1(symbuf);
    size_t nmoved   = 0;
    size_t nr_avail = 0;
    while (tmp_buf1->next != NULL)
    {
        tmp_buf1 = tmp_buf1->next;
        nbuf++;
        nr_avail += tmp_buf1->bufsize;
    }
    // return if all buffers are completely filled
    if (nr_avail == nr)
    {
        printf("t_symtab::clean_symtab(): returning because all buffers are completely filled and hence unfragmented\n");
        return;
    }

    bool fragmented = true;
    while (fragmented)
    {
        fragmented = false;

        // find the last occupied position
        size_t                    tmp_index1   = 0;
        bool                      got_occupied = false;
        std::shared_ptr<t_symbuf> last_occupied(symbuf);
        std::shared_ptr<t_symbuf> ibuf(symbuf);
        while (ibuf->next != NULL)
        {
            ibuf = ibuf->next;
            for (size_t i = 0; i < ibuf->bufsize; i++)
            {
                if (!ibuf->buf[i].expired())
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
        size_t                    tmp_index2 = 0;
        std::shared_ptr<t_symbuf> first_free(symbuf);
        ibuf = symbuf;
        while (ibuf->next != NULL)
        {
            ibuf = ibuf->next;
            for (size_t i = 0; i < ibuf->bufsize; i++)
            {
                // stop at the previously determined last occupied position
                if (i == tmp_index1 && ibuf == last_occupied)
                {
                    break;
                }

                if (ibuf->buf[i].expired())
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
            first_free->buf[tmp_index2].swap(last_occupied->buf[tmp_index1]);
            nmoved++;
        }

    }   // ... end while (fragmented) {}

    // delete empty buffers beginning at the end of the liked list
    tmp_buf1 = symbuf;
    size_t nbufdeleted = 0;
    while (tmp_buf1->next != NULL)
    {
        tmp_buf1 = tmp_buf1->next;
    }
    while (tmp_buf1 != symbuf)
    {
        if (tmp_buf1->buf == NULL)
        {
            break;
        }
        else if (tmp_buf1->buf[0].expired())
        {
            break;
        }
        tmp_buf1 = tmp_buf1->prev.lock();
        tmp_buf1->next.reset();
        nbufdeleted++;
    }

    this->mtx.unlock();

    // update current memory occupied, number of symbol buffers
    // and number of std::string entries
    update_stats();

    // instead of commenting out to avoid compiler warnings because of unused variables
    if (1 == 0)
    {
        printf("%s::clean_symtab():\n", typeid(*this).name());
        printf("moved %li entries, and deleted %li buffers\n", nmoved, nbufdeleted);
        printf("freed %li bytes of memory.\n", tmp_mem - mem);
    }
}

#else

//=====================================================================================
//
// implementations of member functions for the C++98 variant
//
//=====================================================================================

/*! function compares two std::strings, not considering wildcards, and
    returns true if the std::strings are equal and false otherwise
    leading and trailing whitespace is skipped                    */
bool t_symbuf::match_strings(const char * const s1, const char * const s2) const
{
    if (s1 == NULL || s2 == NULL)
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
    while (start1 != NULL && (*start1 == ' ' || *start1 == '\t'))
    {
        start1++;
    }
    while (start2 != NULL && (*start2 == ' ' || *start2 == '\t'))
    {
        start2++;
    }

    // find std::string ends
    const char *end1 = start1;
    const char *end2 = start2;
    // search the end of the std::string marked by the terminating '\0'
    while (end1 != NULL && *end1 != '\0')
    {
        end1++;
    }
    while (end2 != NULL && *end2 != '\0')
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

bool t_symbuf::insert_entry(boost::shared_ptr<boost::shared_ptr<t_symbol> > in)
{
    // Search for a free pointer in the array buf.
    // If found, allocate memory for the new std::string entry
    // and store it in the newly allocated space
    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i] == NULL)
        {
            buf[i] = in;
            return true;
        }
    }

    return false;
}

boost::shared_ptr<boost::shared_ptr<t_symbol> > t_symbuf::get_entry(const char * const s) const
{
    if (s == NULL)
    {
        return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
    }

    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i] == NULL)
        {
            continue;
        }
        else
        {
            boost::shared_ptr<boost::shared_ptr<t_symbol> > tmp_ptr(buf[i]);
            if (match_strings((**tmp_ptr).get_ptr(), s))
            {
                return tmp_ptr;
            }
        }
    }

    return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
}

void t_symbuf::get_stats(size_t * const nr, size_t * const mem)
{
    (*nr)  = 0;
    (*mem) = sizeof(t_symbuf) + sizeof(boost::shared_ptr<boost::shared_ptr<t_symbol> >) * bufsize;
    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i] != NULL)
        {
            (*nr)++;
            size_t      len = std::strlen((**boost::shared_ptr<boost::shared_ptr<t_symbol> >(buf[i])).get_ptr()) + 1;
            long double one_over_secondary_ref_count = 1 /
                static_cast<long double>((**buf[i]).get_ref());
            // size estimate, not guaranteed to be totally accurate
            // because of the shared_ptr contribution
            (*mem) += static_cast<size_t>(
                    (static_cast<long double>(len * sizeof(char) + sizeof(t_symbol))
                     * one_over_secondary_ref_count)
                    + static_cast<long double>(sizeof(boost::shared_ptr<boost::shared_ptr<t_symbol> >)) );
        }
    }
}

void t_symbuf::get_topology(std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> > &topology) const
{

    // topology: std::string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
    typedef boost::shared_ptr<t_symbol>* optr;

    for (size_t i = 0; i < bufsize; ++i)
    {
        if (buf[i] != NULL)
        {
            boost::shared_ptr<boost::shared_ptr<t_symbol> > tmp_ptr(buf[i]);
            std::string this_string((**tmp_ptr).get_ptr());
            // look for an existing entry for this_string
            std::map<std::string, std::map<optr, size_t> >::iterator it1 = topology.find(this_string);
            if (it1 == topology.end())
            {
                // entry not yet in std::map -> create a new entry
                std::map<optr, size_t> tmp_map;
                tmp_map.insert(std::map<optr, size_t>::value_type(tmp_ptr.get(), (tmp_ptr.use_count() - 1)));
                topology.insert(std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> >::value_type(this_string, tmp_map));
            }
            else
            {
                // the entry already exists -> expand it
                // look for an existing entry for tmp_ptr
                std::map<optr, size_t>::iterator it2 = it1->second.find(tmp_ptr.get());
                if (it2 != it1->second.end())
                {
                    // shared ptr set is already present, but the reference count may have changed in the meantime
                    it2->second = tmp_ptr.use_count();
                }
                else
                {
                    // this shared_ptr set is not yet registered, create a new entry
                    it1->second.insert(std::map<optr, size_t>::value_type(tmp_ptr.get(), (tmp_ptr.use_count() - 1)));
                }
            }
        }
    }

}


////// // t_symtab member functions

/*! \brief
    Functions for allocation/deallocation of the symbol table entries,
    which are referenced by shared pointers. In absence of the C++11
    weak_pointer, they also need to make sure that symbol table entries
    are removed from the symbol table as soon as no etzernal reference is left.
    The custom deleters can not be member functions.
 */

//! allocate a new t_symbol on the symtab_inner shared_ptr
t_symbol* symtab_inner_allocator(const size_t &len)
{
    t_symbol* p = new t_symbol;
    p->set_ref(1);
    /* p->set_ptr(new char[len+1]); */
    char *cp = (new char[len+1]);
    cp[0] = '\0';
    p->set_ptr(cp);
    return p;
}

//! deallocate the t_symbol allocated by symtab_inner_allocator
void symtab_inner_deleter(t_symbol* &p)
{
    if (p != NULL)
    {
        if (p->get_ref() == 0)
        {
            const char* tmp_p = p->exchange_ptr(NULL);
            if (tmp_p != NULL)
            {
                delete [] tmp_p;
            }
            delete p;
        }
    }
}

//! allocate a new t_symbol on the symtab_inner shared_ptr
boost::shared_ptr<t_symbol>*  symtab_outer_allocator(const size_t &len)
{
//    void (*del)(t_symbol*&) = symtab_inner_deleter;
    boost::shared_ptr<t_symbol>* p =
        new boost::shared_ptr<t_symbol>( symtab_inner_allocator(len), symtab_inner_deleter );
    return p;
}

//! deallocate the shared_ptr<t_symbol> allocated by symtab_outer_allocator
void symtab_outer_deleter(boost::shared_ptr<t_symbol>* &p)
{
    if (p != NULL)
    {
        // decrement the secondary reference counter,
        // if the resulting secondary reference count
        // is zero, no other symtab_outer pointer referencing
        // the object is present and we can delete it
        if ((*p) != NULL)
        {
            (*p)->decr_ref();
            /*
               size_t rcnt = (*p)->decr_ref();
               if (rcnt == 0)
               {
               }
             */
        }
        delete p;
    }
}


//! trim white space at the std::string ends and remaining chars at the end exceeding maxlen
std::string t_symtab::trim_string(const char * const s1, const size_t maxlen) const
{
    if (s1 == NULL)
    {
        return std::string("");
    }

    int  i = 0,  ii = 0;

    // skip leading whitespace
    while (&s1[i] && (s1[i] == ' ' || s1[i] == '\t'))
    {
        i++;
    }
    // find std::string end (position of the last non-terminating characters)
    ii = i;
    while (&s1[ii+1] != NULL && s1[ii+1] != '\0')
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
    // copy the trimmed std::string to a char array
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

boost::shared_ptr<boost::shared_ptr<t_symbol> > t_symtab::construct_entry(const char * const in, const size_t maxlen, size_t *tmp_mem)
/*!
 * utility function for symbol table, trims whitespace at the beginning of
 * a std::string and trims the std::string to a total maximum length given by maxlen
 * the resulting std::string is written to an output std::string that is managed by a shared_ptr
 * with the correct size
 */
{
    if (in == NULL)
    {
        return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
    }

    size_t      len, i;
    std::string tmp_str(trim_string(in, maxlen));
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
        fprintf(stderr, "Warning in  %s::construct_entry()\n", typeid(*this).name());
        fprintf(stderr, "Truncated exceedingly long std::string to max. length %li. std::string: %s\n", maxlen, in);
        fprintf(stderr, "Error: Insufficient memory for allocating new symbol table buffer!\n");
    }

    void (*odel)(boost::shared_ptr<t_symbol>* &) = symtab_outer_deleter;
    boost::shared_ptr<boost::shared_ptr<t_symbol> > result( symtab_outer_allocator(len), odel );

    char * tmp_cp = (*result)->get_ptr();
    // store the chars forming the std::string
    for (i = 0; i < len; i++)
    {
        tmp_cp[i] = *(s++);
    }
    tmp_cp[i] = '\0';
    // memory occupancy
    (*tmp_mem) = (len + 1) * sizeof(char) + sizeof(t_symbol);
    return result;
}

// the function does not lock this->mtx by itself because it is called
// by other functions that should already own the lock
boost::shared_ptr<boost::shared_ptr<t_symbol> > t_symtab::nonconst_get_entry(const char *s) const
{
    if (s == NULL)
    {
        return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
    }

    boost::shared_ptr<t_symbuf>                     search_pt = symbuf;
    boost::shared_ptr<boost::shared_ptr<t_symbol> > result;

    while (search_pt != NULL)
    {
        result = search_pt->get_entry(s);
        if (result != NULL)
        {
            return result;
        }
        search_pt = search_pt->next;
    }

    return result;
}

boost::shared_ptr<boost::shared_ptr<t_symbol> > t_symtab::nonconst_insert_entry(const char * const s)
{
    if (s == NULL)
    {
        return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
    }

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    gmx::lock_guard<gmx::Mutex> lock(this->mtx);

    // first try if the entry exists already, otherwise create a new entry
    boost::shared_ptr<boost::shared_ptr<t_symbol> > result(nonconst_get_entry(s));
    if (result != NULL)
    {
        return result;
    }
    else
    {

        // look for an existing buffer with an unoccupied entry
        // and if it exists store the new entry there
        boost::shared_ptr<t_symbuf> this_pt = symbuf;
        boost::shared_ptr<t_symbuf> prev_pt = boost::shared_ptr<t_symbuf>();
        size_t                      tmp_mem;

        result = construct_entry(s, max_entry_size, &tmp_mem);
        bool success = false;
        while (this_pt != NULL)
        {
            try
            {
                success = this_pt->insert_entry(result);
            }
            catch (std::logic_error &)
            {
                fprintf(stderr, "[exception caught] in %s::insert_entry():\n", typeid(*this).name());
                fprintf(stderr, "failed to insert entry because of internal inconsistency\n");
                throw;
            }
            catch (std::bad_alloc &)
            {
                fprintf(stderr, "[exception caught] in %s::insert_entry():\n", typeid(*this).name());
                fprintf(stderr, "failed to allocate memory for std::string entry\n");
                throw;
            }

            if (success)
            {
                nr++;
                mem += tmp_mem;
                return result;
            }
            prev_pt = this_pt;
            this_pt = this_pt->next;
        }

        // allocate memory for a new buffer because
        // all preexisting buffers were full
        if (this_pt == NULL)
        {
            this_pt = prev_pt;
            // we could still have no t_symbuf at all
            if (this_pt == NULL)
            {
                symbuf  = boost::shared_ptr<t_symbuf>(new t_symbuf());
                this_pt = symbuf;
            }
            else
            {
                this_pt->next = boost::shared_ptr<t_symbuf>(new t_symbuf());
                this_pt       = this_pt->next;
            }
            this_pt->next.reset();
            this_pt->prev = prev_pt;

            try
            {
                this_pt->bufsize = bufsize;
                this_pt->alloc_symbuf(&tmp_mem);
            }
            catch (std::bad_alloc &)
            {
                fprintf(stderr, "Error in  %s::insert_entry()\n", typeid(*this).name());
                fprintf(stderr, "Insufficient memory for allocating new symbol table buffer!\n");
                throw;
            }
            ;
            nbuf++;
            mem += tmp_mem;

        }
        else
        {
            fprintf(stderr, "Error in  %s::insert_entry()\n", typeid(*this).name());
            fprintf(stderr, "Failed to find end of linked list for allocating a new symbol table buffer for unknown reason.\n");
            throw std::logic_error("Internal error in t_symtab");
        }

        // store the new std::string entry
        success = this_pt->insert_entry(result);

        nr++;
        mem += tmp_mem;

        // check whether the entry was inserted successfully
        if (!success)
        {
            // this should never happen as long as the code of t_symtab and
            // t_symbuf is internally consistent
            fprintf(stderr, "Error in  %s::insert_entry()\n", typeid(*this).name());
            fprintf(stderr, "Could not insert new entry for unknown reason.\n");
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
        return result;
    }

//// if (result) {
////   std::cout << "  1st level is OK. use count " << result.use_count() << std::endl;
////   if (*result) {
////     std::cout << "    2nd level is OK. use count " << result.use_count() << std::endl;
////     std::cout << "    (*result)->get_ref()       : " << (*result)->get_ref() << std::endl;
////     std::cout << "    (*result)->get_ptr()       : " << static_cast<void*>((*result)->get_ptr()) << std::endl;
////     if ((*result)->get_ptr()) {
////       std::cout << "    (*result)->get_ptr() std::string: " << (*result)->get_ptr() << std::endl;
////     }
////   } else {
////     std::cout << "    2nd level is NOT OK." << std::endl;
////   }
//// } else {
////   std::cout << "  1st level is NOT OK." << std::endl;
//// }

    return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
}

boost::shared_ptr<boost::shared_ptr<t_symbol> > t_symtab::nonconst_rename_entry(const char * const before, const char * const after)
{
    if (before == NULL || after == NULL)
    {
        return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
    }

    // the nesting of multiple member functions necessitates manual locking
    // and unlocking instead of the convenient use of lock_guard
    // changes by other threads could take place in between the periods of
    // owning the lock but that should not cause problems as long as there
    // are no simultaneous calls to the function by multiple threads
    // first try whether the entry before exists already,
    // otherwise delegate the job to insert_entry()
    this->mtx.lock();
    boost::shared_ptr<boost::shared_ptr<t_symbol> > result(nonconst_get_entry(before));
    this->mtx.unlock();

    if (result)
    {

        // If the entry "after" already exists, exchange the target of the
        // pointer pointing to "before". Otherwise, store the new std::string
        // "after" in a newly allocated memory area and let the pointer
        // to "before" point to it.
        boost::shared_ptr<boost::shared_ptr<t_symbol> > tmp_ptr(nonconst_insert_entry(after));
        if (tmp_ptr)
        {

            // other threads might have made changes that removed the
            // previously found entry for before in the meantime
            // "before" is preserved by the referencing shared_ptr instance "result"
            // "after", by "tmp_ptr"
            gmx::lock_guard<gmx::Mutex> lock(this->mtx);
            result = nonconst_get_entry(before);
            while (result)
            {
                // There seems to be no way to manipulate the control block of the shared pointers
                // such that two sets of shared_pointer instances can be merged, a custom shared
                // pointer with that capability should in principle be possible but detach the code
                // from possible future improvements of std::shared_ptr. Is there a safe, strictly
                // portable way via inheritance from shared_ptr? probably not.
                (*tmp_ptr)->incr_ref();
                boost::shared_ptr<t_symbol> save_last(*result);
                *result = *tmp_ptr;
                result  = nonconst_get_entry(before);
                if (!result)
                {
                    save_last->set_ref((size_t)0);
                }
            }

            result = nonconst_get_entry(after);
            return result;
        }
        else
        {
            // this should never happen as long as the code of t_symtab and
            // t_symbuf is internally consistent
            fprintf(stderr, "Error in  %s::nonconst_rename_entry()\n", typeid(*this).name());
            fprintf(stderr, "Could not insert or locate entry for unknown reason.\n");
            throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
        }
    }
    else
    {
        // not a fatal error, simply create a new entry or return a shared_ptr
        // to a possible preexisting entry for "after"
        return nonconst_insert_entry(after);
    }

    return boost::shared_ptr<boost::shared_ptr<t_symbol> >();
}

void t_symtab::get_string_stats(const char * const s, size_t *ninst, size_t *nptr) const
{
    if (s == NULL)
    {
        return;
    }

    // topology: std::string -> unique addresses of ptrs managed by distinct shared pointer sets
    //                   -> reference counts of the shared pointer sets
    typedef std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> > map_type;
    typedef std::map<boost::shared_ptr<t_symbol>*, size_t> imap_type;
    map_type           topology;
    get_topology(topology);
    map_type::iterator it1 = topology.find(std::string(s));

    *nptr  = 0;
    *ninst = 0;
    if (it1 != topology.end())
    {

        for (imap_type::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
            (*nptr)++;
            //! one less because the shared_ptr instance in the symbol table also counts
            (*ninst) += it2->second - 1;
        }

        if (*nptr != (*nonconst_get_entry(s))->get_ref())
        {
            // this should never happen as long as the code of t_symtab and
            // t_symbuf is internally consistent
            fprintf(stderr, "Error in  %s::get_string_stats()\n", typeid(*this).name());
            fprintf(stderr, "Detected inconsistent secondary reference counts.\n");
            fprintf(stderr, " get_topology(): %li vs. t_symbol: %li\n", *nptr, (*nonconst_get_entry(s))->get_ref());
        }

    }
    else
    {
        // this should never happen as long as the code of t_symtab and
        // t_symbuf is internally consistent
        fprintf(stderr, "Error in  %s::get_string_stats()\n", typeid(*this).name());
        fprintf(stderr, "Could not find entry: \"%s\"\n", s);
        throw std::logic_error("Internal error in " + std::string(typeid(*this).name()));
    }

}

void t_symtab::get_topology(std::map<std::string, std::map<boost::shared_ptr<t_symbol>*, size_t> > &topology) const
{

    boost::shared_ptr<t_symbuf> search_pt = symbuf;

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    gmx::lock_guard<gmx::Mutex> lock(this->mtx);

    while (search_pt != NULL)
    {
        search_pt->get_topology(topology);
        search_pt = search_pt->next;
    }

}

void t_symtab::update_stats()
{
    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    gmx::lock_guard<gmx::Mutex> lock(this->mtx);
    mem  = sizeof(t_symtab);
    nr   = 0;
    nbuf = 0;

    // only proceed if there is any symbol buffer yet
    if (symbuf == NULL)
    {
        return;
    }

    // sum up the numbers of t_symbufs and their std::string entries
    // and the memory occupied by them
    size_t tmp_mem = 0;
    size_t tmp_nr  = 0;
    boost::shared_ptr<t_symbuf> tmp_buf1(symbuf);
    while (tmp_buf1 != NULL)
    {
        tmp_buf1->get_stats(&tmp_nr, &tmp_mem);
        nr  += tmp_nr;
        mem += tmp_mem;
        nbuf++;
        tmp_buf1 = tmp_buf1->next;
    }
}

void t_symtab::clean()
{
    // return if there is no stored data at all
    if (symbuf == NULL)
    {
        printf("t_symtab::clean_symtab(): returning because there is no data allocated\n");
        return;
    }

    // determine number of entries and symbol buffers
    update_stats();
    size_t tmp_mem = mem;

    // the mutex is locked until the end of the scope
    // and automatically released when leaving the scope
    this->mtx.lock();

    // remove all entries that are not referenced anymore from outside
    // C++11 weak_pointers in t_symtab would solve that automatically)
    if (symbuf != NULL)
    {
        boost::shared_ptr<t_symbuf> tmp_buf1(symbuf);
        while (tmp_buf1 != NULL)
        {
            tmp_buf1->remove_orphan_entries();
            tmp_buf1 = tmp_buf1->next;
        }
    }

    boost::shared_ptr<t_symbuf>   tmp_buf1(symbuf);
    size_t nmoved      = 0;
    size_t nr_avail    = 0;
    while (tmp_buf1->next != NULL)
    {
        tmp_buf1 = tmp_buf1->next;
        nbuf++;
        nr_avail += tmp_buf1->bufsize;
    }
    // return if all buffers are completely filled
    if (nr_avail == nr)
    {
        printf("t_symtab::clean_symtab(): returning because all buffers are completely filled and hence unfragmented\n");
        return;
    }

    bool fragmented = true;
    while (fragmented)
    {
        fragmented = false;

        // find the last occupied position
        size_t tmp_index1   = 0;
        bool   got_occupied = false;
        boost::shared_ptr<t_symbuf> last_occupied(symbuf);
        boost::shared_ptr<t_symbuf> ibuf(symbuf);
        while (ibuf->next != NULL)
        {
            ibuf = ibuf->next;
            for (size_t i = 0; i < ibuf->bufsize; i++)
            {
                if (ibuf->buf[i] != NULL)
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
        size_t tmp_index2 = 0;
        boost::shared_ptr<t_symbuf> first_free(symbuf);
        ibuf = symbuf;
        while (ibuf->next != NULL)
        {
            ibuf = ibuf->next;
            for (size_t i = 0; i < ibuf->bufsize; i++)
            {
                // stop at the previously determined last occupied position
                if (i == tmp_index1 && ibuf == last_occupied)
                {
                    break;
                }

                if (ibuf->buf[i] == NULL)
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
            first_free->buf[tmp_index2].swap(last_occupied->buf[tmp_index1]);
            nmoved++;
        }

    }   // ... end while (fragmented) {}

    // delete empty buffers beginning at the end of the liked list
    tmp_buf1 = symbuf;
    size_t nbufdeleted = 0;
    while (tmp_buf1->next != NULL)
    {
        tmp_buf1 = tmp_buf1->next;
    }
    while (tmp_buf1 != symbuf)
    {
        if (tmp_buf1->buf == NULL)
        {
            break;
        }
        else if (tmp_buf1->buf[0] == NULL)
        {
            break;
        }
        tmp_buf1 = tmp_buf1->prev;
        tmp_buf1->next->prev.reset();
        tmp_buf1->next.reset();
        nbufdeleted++;
    }

    this->mtx.unlock();

    // update current memory occupied, number of symbol buffers
    // and number of std::string entries
    update_stats();

    // instead of commenting out to avoid compiler warnings because of unused variables
    if (1 == 0)
    {
        printf("%s::clean_symtab():\n", typeid(*this).name());
        printf("moved %li entries, and deleted %li buffers\n", nmoved, nbufdeleted);
        printf("freed %li bytes of memory.\n", tmp_mem - mem);
    }
}

#endif

} // namespace gmx
