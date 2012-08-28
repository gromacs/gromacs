#pragma once

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
    msg_impl(MsgTy& v, int tag = 0) : m_data(v), m_tag(tag){ }

    // Returns the address to the first element of the contained data
    inline void* addr() {
        return (void*) mpi_type_traits<MsgTy>::get_addr(m_data);
    }
    inline const void* addr() const {
        return (const void*) mpi_type_traits<MsgTy>::get_addr(m_data);
    }

    inline MsgTy& get() { return m_data; }
    inline const MsgTy& get() const { return m_data; }

    // Returns the dimension of this message
    inline size_t size() const {
        return mpi_type_traits<MsgTy>::get_size(m_data);
    }

    inline MPI_Datatype type() const {
        return mpi_type_traits<MsgTy>::get_type(m_data);
    }

    inline MPI_Datatype type_range(int begin, int end) const {
        return mpi_type_traits<MsgTy>::get_type_range(m_data, begin, end);
    }

    inline MPI_Datatype type_selection (const std::vector<int>& sel, int begin, int end) const {
        return mpi_type_traits<MsgTy>::get_type_selection(m_data, sel, begin, end);
    }

    inline bool type_static() const {
        return mpi_type_traits<MsgTy>::get_static();
    }

    // getter/setter for m_tag
    inline const int& tag() const { return m_tag; }
    inline int& tag() { return m_tag; }

private:
    MsgTy&  m_data;
    int     m_tag;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Specializaton for class MsgTy for const types in this case we don't keep the
// reference to the object passed to the constructor, but we make a copy of it
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class MsgTy>
struct msg_impl <const MsgTy> {
    typedef const MsgTy value_type;

    // Builds a msg wrapping v
    msg_impl(const MsgTy& v, int tag = 0) : m_data(v), m_tag(tag){ }

    // Returns the enclosed data
    inline const void* addr() const {
        return mpi_type_traits<MsgTy>::get_addr(m_data);
    }

    inline const MsgTy& get() const { return m_data; }

    // Returns the dimension of this message
    inline size_t size() const {
        return mpi_type_traits<MsgTy>::get_size(m_data);
    }

    inline MPI_Datatype type() const {
        return mpi_type_traits<MsgTy>::get_type(m_data);
    }

    // getter/setter for m_tag
    inline const int& tag() const { return m_tag; }
    inline int& tag() { return m_tag; }

private:
    const MsgTy m_data;
    int         m_tag;
};

template <class T>
msg_impl<T> msg(T& raw, int tag=0) { return msg_impl<T>(raw, tag); }
}
