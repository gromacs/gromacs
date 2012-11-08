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
 *
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
/*! \file
 * \brief
 * Declares derived iterators identical to those in boost.
 *
 * \author Roland Schulz <roland@utk.edu>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ITERATOR_H
#define GMX_UTILITY_ITERATOR_H
#ifdef USE_BOOST_ITERATORS //can be manually set to test compatibility
#include <boost/iterator/permutation_iterator.hpp>
namespace gmx {
using boost::permutation_iterator;
using boost::make_permutation_iterator;
}
#else
#include <algorithm>
#include <iterator>
namespace gmx {
template< class ElementIterator, class IndexIterator>
class permutation_iterator : public std::iterator<
        typename IndexIterator::iterator_category,
        typename ElementIterator::value_type,
        typename IndexIterator::difference_type,
        typename ElementIterator::pointer,
        typename ElementIterator::reference>
{
    //TODO: is there some way to not have to list all types twice?
    typedef std::iterator<
        typename IndexIterator::iterator_category,
        typename ElementIterator::value_type,
        typename IndexIterator::difference_type,
        typename ElementIterator::pointer,
        typename ElementIterator::reference> BaseIterator;
public:
    permutation_iterator() : m_elt_iter(), m_idx_iter()  {}

    permutation_iterator(ElementIterator x, IndexIterator y)
            : m_elt_iter(x), m_idx_iter(y) {}

    //returns *(m_elt_iter+*m_idx_iter)
    typename BaseIterator::reference operator*()
    {
        ElementIterator curr = m_elt_iter;
        std::advance(curr, *m_idx_iter);
        return *curr;
    }
    permutation_iterator& operator++() { ++m_idx_iter; return *this; }
    permutation_iterator& operator+=(typename BaseIterator::difference_type i)
    {
        m_idx_iter+=i; return *this;
    }
    typename BaseIterator::difference_type operator-(
            const permutation_iterator &other)
    {
        return m_idx_iter-other.m_idx_iter;
    }
    bool operator==(const permutation_iterator &other)
    {
        assert(m_elt_iter==other.m_elt_iter);
        return m_idx_iter==other.m_idx_iter;
    }
    bool operator!=(const permutation_iterator &other)
    {
        return !(*this==other);
    }
private:
    ElementIterator m_elt_iter;
    IndexIterator m_idx_iter;
};
//! Creates a permutation iterator
template <class ElementIterator, class IndexIterator>
permutation_iterator<ElementIterator, IndexIterator>
make_permutation_iterator( ElementIterator e, IndexIterator i )
{
    return permutation_iterator<ElementIterator, IndexIterator>( e, i );
}
} //namespace gmx
#endif
#endif
