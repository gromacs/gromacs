/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */

#pragma once

#include "decls.h"

#include "type_traits.h"

namespace mpi
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// msg: represent a single message which can be provided to the <<, <, >>, >
// operations
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class MsgTy>
struct msg_impl {
    typedef MsgTy value_type;

    // Builds a msg wrapping v
    msg_impl(value_type& v, int tag = 0) : m_data(v), m_tag(tag) {
        m_dt = get_type_layout(v);
    }

    // Returns the address to the first element of the contained data
    inline void* addr() {
        return m_dt.get_ptr();
    }
    inline const void* addr() const {
        return m_dt.get_ptr();
    }

    inline int size() const {
        return m_dt.get_size();
    }

    inline MsgTy& get() { return m_data; }
    inline const MsgTy& get() const { return m_data; }

    inline MPI_Datatype type() const {
        return m_dt.get();
    }

    // getter/setter for m_tag
    inline const int& tag() const { return m_tag; }
    inline int& tag() { return m_tag; }

private:
    MsgTy&   m_data;
    data_layout m_dt;
    int      m_tag;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Specializaton for class MsgTy for const types in this case we don't keep the
// reference to the object passed to the constructor, but we make a copy of it
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class MsgTy>
struct msg_impl <const MsgTy> {
    typedef const MsgTy value_type;

    // Builds a msg wrapping v
    msg_impl(const MsgTy& v, int tag = 0) : m_data(v), m_tag(tag){
        m_dt = get_type_layout(v);
    }

    // Returns the enclosed data
    inline const void* addr() const {
        return m_dt.get_ptr();
    }

    inline int size() const {
        return m_dt.get_size();
    }

    inline const MsgTy& get() const { return m_data; }

    inline MPI_Datatype type() const {
        return m_dt.get();
    }

    // getter/setter for m_tag
    inline const int& tag() const { return m_tag; }
    inline int& tag() { return m_tag; }

private:
    const MsgTy m_data;
    data_layout m_dt;
    int         m_tag;
};

template <class T>
msg_impl<T> msg(T& raw, int tag=0) { return msg_impl<T>(raw, tag); }

}
