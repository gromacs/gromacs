/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
// Copyright 2017-2019 Martin Moene
//
// https://github.com/martinmoene/string-view-lite
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
/*! \internal \file
 * \brief Tests for C++14-compatible implementation of std::string_view.
 *
 * These tests are similer to those used in commit bf5824916b6895ccab0dbc2431520ee3b6d4f27f of
 * https://github.com/martinmoene/string_view-lite.git. The code has not
 * been linted with uncrustify so that any future updates to this
 * active repo can be incorporated easily, and //NOLINT comments added
 * to suppress clang-tidy warnings. The form of those
 * changes have been made to simplify the contents, while making it
 * easy to import any bug fixes that may appear in the source
 * repository.
 *
 * \todo Remove when requiring C++17, which has a standardized version
 * of std::string_view.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_compat
 */
#include "gmxpre.h"

#include "gromacs/compat/string_view.h"

#include <initializer_list>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/strconvert.h"

// GMX modification to suppress Doxygen checking
#ifndef DOXYGEN

namespace gmx
{

using namespace compat;

namespace {

template< class T >
T * data( std::vector<T> & v )
{
#if nssv_CPP11_OR_GREATER
    return v.data();
#else
    return &v[0];
#endif
}

#define nssv_STD_SV_OR( expr )  ( nssv_USES_STD_STRING_VIEW || (expr) )

typedef string_view::size_type size_type;

// 24.4.2.1 Construction and assignment:

TEST(StringViewTest, CanDefaultConstructEmptyView)
{
    SCOPED_TRACE( "string_view: Allows to default construct an empty string_view" );

    string_view sv;

    // use parenthesis with data() to prevent lest from using sv.data() as C-string:

    EXPECT_TRUE(  sv.size() == size_type( 0 ) );
    EXPECT_TRUE( (sv.data() == nssv_nullptr)  );
}

TEST(StringViewTest, CanConstructFromPointerAndSize)
{
    SCOPED_TRACE( "string_view: Allows to construct from pointer and size" );

    string_view sv( "hello world", 5 );

    EXPECT_TRUE(   sv.size() == size_type( 5 ) );
    EXPECT_TRUE( *(sv.data() + 0) == 'h'       );
    EXPECT_TRUE( *(sv.data() + 4) == 'o'       );
}

TEST(StringViewTest, CanConstructFromCString)
{
    SCOPED_TRACE( "string_view: Allows to construct from C-string" );

    string_view sv( "hello world" );

    EXPECT_TRUE(   sv.size() == size_type( 11 ) );
    EXPECT_TRUE( *(sv.data() +  0) == 'h'       );
    EXPECT_TRUE( *(sv.data() + 10) == 'd'       );
}

TEST(StringViewTest, CanCopyConstructFromEmptyView)
{
    SCOPED_TRACE( "string_view: Allows to copy-construct from empty string_view" );

    string_view sv1;

    string_view sv2( sv1 );

    // use parenthesis with data() to prevent lest from using sv.data() as C-string:

    EXPECT_TRUE(  sv2.size() == size_type( 0 ) );
    EXPECT_TRUE( (sv2.data() == nssv_nullptr)  );
}

TEST(StringViewTest, CanCopyConstructFromNonEmptyView)
{
    SCOPED_TRACE( "string_view: Allows to copy-construct from non-empty string_view" );

    string_view sv1( "hello world", 5 );

    string_view sv2( sv1 );

    EXPECT_TRUE(   sv2.size()      == sv1.size()  );
    EXPECT_TRUE(  (sv2.data()      == sv1.data()) );
    EXPECT_TRUE( *(sv2.data() + 0) == 'h'         );
    EXPECT_TRUE( *(sv2.data() + 4) == 'o'         );
}

// Assignment:

TEST(StringViewTest, CanCopyAssingFromEmptyView)
{
    SCOPED_TRACE( "string_view: Allows to copy-assign from empty string_view" );

    string_view sv1;
    string_view sv2;

    sv2 = sv1;

    // use parenthesis with data() to prevent lest from using sv.data() as C-string:

    EXPECT_TRUE(  sv2.size() == size_type( 0 ) );
    EXPECT_TRUE( (sv2.data() == nssv_nullptr) );
}

TEST(StringViewTest, CanCopyAssingFromNonEmptyView)
{
    SCOPED_TRACE( "string_view: Allows to copy-assign from non-empty string_view" );

    string_view sv1( "hello world", 5 );
    string_view sv2;

    sv2 = sv1;

    // use parenthesis with data() to prevent lest from using sv.data() as C-string:

    EXPECT_TRUE(   sv2.size()      == sv1.size()  );
    EXPECT_TRUE(  (sv2.data()      == sv1.data()) );
    EXPECT_TRUE( *(sv2.data() + 0) == 'h'         );
    EXPECT_TRUE( *(sv2.data() + 4) == 'o'         );
}

// 24.4.2.2 Iterator support:

TEST(StringViewTest, AllowForwardIteration)
{
    SCOPED_TRACE( "string_view: Allows forward iteration" );

    char hello[] = "hello";
    string_view sv( hello );

    for ( string_view::iterator pos = sv.begin(); pos != sv.end(); ++pos )
    {
        typedef std::iterator_traits< string_view::iterator >::difference_type difference_type;

        difference_type i = std::distance( sv.begin(), pos );
        EXPECT_TRUE( *pos == *(sv.data() + i) );
    }
}

TEST(StringViewTest, AllowConstForwardIteration)
{
    SCOPED_TRACE( "string_view: Allows const forward iteration" );

    char hello[] = "hello";
    string_view sv( hello );

    for ( string_view::const_iterator pos = sv.begin(); pos != sv.end(); ++pos )
    {
        typedef std::iterator_traits< string_view::const_iterator >::difference_type difference_type;

        difference_type i = std::distance( sv.cbegin(), pos );
        EXPECT_TRUE( *pos == *(sv.data() + i) );
    }
}

TEST(StringViewTest, AllowReverseIteration)
{
    SCOPED_TRACE( "string_view: Allows reverse iteration" );

    char hello[] = "hello";
    string_view sv( hello );

    for ( string_view::reverse_iterator pos = sv.rbegin(); pos != sv.rend(); ++pos )
    {
        typedef std::iterator_traits< string_view::reverse_iterator >::difference_type difference_type;

        difference_type dist = std::distance( sv.rbegin(), pos );
        EXPECT_TRUE( *pos == *(sv.data() + sv.size() - 1 - dist) );
    }
}

TEST(StringViewTest, AllowConstReverseIteration)
{
    SCOPED_TRACE( "string_view: Allows const reverse iteration" );

    char hello[] = "hello";
    string_view sv( hello );

    for ( string_view::const_reverse_iterator pos = sv.crbegin(); pos != sv.crend(); ++pos )
    {
        typedef std::iterator_traits< string_view::const_reverse_iterator >::difference_type difference_type;

        difference_type  dist = std::distance( sv.crbegin(), pos );
        EXPECT_TRUE( *pos == *(sv.data() + sv.size() - 1 - dist) );
    }
}

// 24.4.2.3 Capacity:

TEST(StringViewTest, CanObtainSizeFromViewViaSize)
{
    SCOPED_TRACE( "string_view: Allows to obtain the size of the view via size()" );

    char hello[] = "hello";
    string_view sv( hello );

    EXPECT_TRUE( sv.size() == std::char_traits<char>::length(hello) );
}

TEST(StringViewTest, CanObtainSizeFromViewViaLength)
{
    SCOPED_TRACE( "string_view: Allows to obtain the size of the view via length()" );

    char hello[] = "hello";
    string_view sv( hello );

    EXPECT_TRUE( sv.length() == std::char_traits<char>::length(hello) );
}

TEST(StringViewTest, CanObtainMaxSizeViaMaxSize)
{
    SCOPED_TRACE( "string_view: Allows to obtain the maximum size a view can be via max_size()" );

    // "large"
    EXPECT_TRUE( string_view().max_size() >= (std::numeric_limits< string_view::size_type >::max)() / 10 );
}

TEST(StringViewTest, CanCheckForEmptyStringWithEmpty)
{
    SCOPED_TRACE( "string_view: Allows to check for an empty string_view via empty()" );

    string_view sve;
    string_view svne("hello");

    EXPECT_TRUE(      sve.size() == size_type( 0 ) );
    EXPECT_TRUE(      sve.empty() );
    EXPECT_FALSE( svne.empty() );
}

// 24.4.2.4 Element access:

TEST(StringViewTest, CanAccessElementViaArrayIndex)
{
    SCOPED_TRACE( "string_view: Allows to observe an element via array indexing" );

    // Requires: index < sv.size()

    char hello[] = "hello";
    string_view sv( hello );

    for ( size_type i = 0; i < sv.size(); ++i )
    {
        EXPECT_TRUE( sv[i] == hello[i] );
    }
}

TEST(StringViewTest, CanAccessElementViaAt)
{
    SCOPED_TRACE( "string_view: Allows to observe an element via at()" );

    char hello[] = "hello";
    string_view sv( hello );

    for ( size_type i = 0; i < sv.size(); ++i )
    {
        EXPECT_TRUE( sv.at(i) == hello[i] );
    }
}

TEST(StringViewTest, ThrowsOnOutOfBoundsAccess)
{
    SCOPED_TRACE( "string_view: Throws at observing an element via at() with an index of size() or larger" );

    string_view sv("hello");

    EXPECT_ANY_THROW(   sv.at( sv.size()     ) );
    EXPECT_NO_THROW( sv.at( sv.size() - 1 ) );
}

TEST(StringViewTest, CanAccessAllElementsViaData)
{
    SCOPED_TRACE( "string_view: Allows to observe elements via data()" );

    char hello[] = "hello";
    string_view sv( hello );

    EXPECT_TRUE( *sv.data() == *sv.begin() );

    for ( size_type i = 0; i < sv.size(); ++i )
    {
        EXPECT_TRUE( sv.data()[i] == hello[i] );
    }
}

TEST(StringViewTest, DataFromEmptyStringIsNullptr)
{
    SCOPED_TRACE( "string_view: Yields nullptr (or NULL) with data() for an empty string_view" );

    string_view sv;

    // use parenthesis with data() to prevent lest from using sv.data() as C-string:

    EXPECT_TRUE( (sv.data() == nssv_nullptr) );
}

// 24.4.2.5 Modifiers:

TEST(StringViewTest, CanRemovePrefix)
{
    SCOPED_TRACE( "string_view: Allows to remove a prefix of n elements" );

    char hello[] = "hello world";
    string_view sv( hello );

    sv.remove_prefix( 6 );

    EXPECT_TRUE( sv.size() == size_type( 5 ) );
    EXPECT_TRUE( std::equal( sv.begin(), sv.end(), hello + 6) );
}

TEST(StringViewTest, CanRemoveSuffix)
{
    SCOPED_TRACE( "string_view: Allows to remove a suffix of n elements" );

    char hello[] = "hello world";
    string_view sv( hello );

    sv.remove_suffix( 6 );

    EXPECT_TRUE( sv.size() == size_type( 5 ) );
    EXPECT_TRUE( std::equal( sv.begin(), sv.end(), hello ) );
}

TEST(StringViewTest, CanSwapWithOtherView)
{
    SCOPED_TRACE( "string_view: Allows to swap with other string_view" );

    char hello[]  = "hello";
    char world[]  = "world";
    string_view sv1( hello );
    string_view sv2( world );

    sv1.swap( sv2 );

    EXPECT_TRUE( std::equal( sv1.begin(), sv1.end(), world )  );
    EXPECT_TRUE( std::equal( sv2.begin(), sv2.end(), hello )  );
}

// 24.4.2.6 String operations:

TEST(StringViewTest, CanCopySubstringWithCopy)
{
    SCOPED_TRACE( "string_view: Allows to copy a substring of length n, starting at position pos (default: 0) via copy()" );

    char hello[]  = "hello world";
    string_view sv( hello );

    {
        std::vector<string_view::value_type> vec( sv.size() );

        sv.copy( data(vec), vec.size() );

        EXPECT_TRUE( std::equal( vec.begin(), vec.end(), hello )  );
    }{
        std::size_t offset = 3U; std::size_t length = 4U;
        std::vector<string_view::value_type> vec( length );

        sv.copy( data(vec), length, offset );

        EXPECT_TRUE( std::equal( vec.begin(), vec.end(), hello + offset )  );
    }
}

TEST(StringViewTest, ThrowsOnOutOfBoundsCopy)
{
    SCOPED_TRACE( "string_view: Throws if requested position of copy() exceeds string_view's size()" );

    string_view sv("hello world");
    std::vector<string_view::value_type> vec( sv.size() );

    EXPECT_ANY_THROW(   sv.copy( data(vec), sv.size() - 4, sv.size() + 1 ) );
    EXPECT_NO_THROW( sv.copy( data(vec), sv.size() - 4, sv.size() + 0 ) );
}

TEST(StringViewTest, CanObtainSubstringWithSubstr)
{
    SCOPED_TRACE( "string_view: Allow to obtain a sub string, starting at position pos (default: 0) and of length n (default full), via substr()" );

    char hello[] = "hello world";
    string_view sv( hello );

    {
#if nssv_USES_STD_STRING_VIEW && defined(__GNUC__)
        EXPECT_TRUE( !!"std::string_view::substr(), with default pos, count is not available with GNU C" );
#else
        EXPECT_TRUE( std::equal( sv.begin(), sv.end(), sv.substr().begin() ) );
#endif
    }{
        string_view subv = sv.substr( 6 );

        EXPECT_TRUE( std::equal( subv.begin(), subv.end(), hello + 6 ) );
    }{
        string_view subv = sv.substr( 3, 4 );

        EXPECT_TRUE( std::equal( subv.begin(), subv.end(), hello + 3 ) );
    }
}

TEST(StringViewTest, ThrowsOnOutOfBoundsSubstr)
{
    SCOPED_TRACE( "string_view: Throws if requested position of substr() exceeds string_view's size()" );

    string_view sv("hello world");

    EXPECT_ANY_THROW(   sv.substr( sv.size() + 1 ) );
    EXPECT_NO_THROW( sv.substr( sv.size() + 0 ) );
}

TEST(StringViewTest, CanLexicallyCompareViewWithCompare)
{
    SCOPED_TRACE( "string_view: Allows to lexically compare to another string_view via compare(), (1)" );

    char hello[] = "hello";
    char world[] = "world";

    EXPECT_TRUE( string_view( hello ).compare( string_view( hello ) ) == 0 );
    EXPECT_TRUE( string_view( hello ).compare( string_view( world ) ) <  0 );
    EXPECT_TRUE( string_view( world ).compare( string_view( hello ) ) >  0 );

    char hello_sp[] = "hello ";

    EXPECT_TRUE( string_view( hello    ).compare( string_view( hello_sp ) ) < 0 );
    EXPECT_TRUE( string_view( hello_sp ).compare( string_view( hello    ) ) > 0 );
}

TEST(StringViewTest, CanCompareEmptyViewsWIthCompare)
{
    SCOPED_TRACE( "string_view: Allows to compare empty string_views as equal via compare(), (1)" );

    EXPECT_TRUE( string_view().compare( string_view() ) == 0 );
}

TEST(StringViewTest, CanCompareSubStringWithViewViaCompare)
{
    SCOPED_TRACE( "string_view: Allows to compare a sub string to another string_view via compare(), (2)" );

    string_view sv1("hello world");
    string_view sv2("world");

    EXPECT_TRUE( sv1.compare ( 0, sv1.length(), sv1 ) == 0 );
    EXPECT_TRUE( sv1.compare ( 6, 5, sv2 ) == 0 );
    EXPECT_TRUE( sv1.compare ( 0, 5, sv2 ) <  0 );
    EXPECT_TRUE( sv2.compare ( 0, 5, sv1 ) >  0 );
}

TEST(StringViewTest, CanCompareSubStringWithSubStringViewViaCompare)
{
    SCOPED_TRACE( "string_view: Allows to compare a sub string to another string_view sub string via compare(), (3)" );

    string_view sv1("hello world");

    EXPECT_TRUE( sv1.compare ( 0, sv1.length(), sv1 ) == 0 );
    EXPECT_TRUE( sv1.compare ( 6, 5, sv1, 6, 5 ) ==  0 );
    EXPECT_TRUE( sv1.compare ( 0, 5, sv1, 6, 5 )  <  0 );
    EXPECT_TRUE( sv1.compare ( 6, 5, sv1, 0, 5 )  >  0 );
}

TEST(StringViewTest, CanCompareToCStringViaCompare)
{
    SCOPED_TRACE( "string_view: Allows to compare to a C-string via compare(), (4)" );

    char hello[] = "hello";
    char world[] = "world";

    EXPECT_TRUE( string_view( hello ).compare( hello ) == 0 );
    EXPECT_TRUE( string_view( hello ).compare( world ) <  0 );
    EXPECT_TRUE( string_view( world ).compare( hello ) >  0 );
}

TEST(StringViewTest, CanCompareSubStringToCStringViaCompare)
{
    SCOPED_TRACE( "string_view: Allows to compare a sub string to a C-string via compare(), (5)" );

    char hello[] = "hello world";
    char world[] = "world";

    EXPECT_TRUE( string_view( hello ).compare( 6, 5, world ) == 0 );
    EXPECT_TRUE( string_view( hello ).compare( world ) <  0 );
    EXPECT_TRUE( string_view( world ).compare( hello ) >  0 );
}

TEST(StringViewTest, CanCompareSubStringToCStringPrefixViaCompare)
{
    SCOPED_TRACE( "string_view: Allows to compare a sub string to a C-string prefix via compare(), (6)" );

    char hello[] = "hello world";
    char world[] = "world";

    EXPECT_TRUE( string_view( hello ).compare( 6, 5, world, 5 ) == 0 );
    EXPECT_TRUE( string_view( hello ).compare( 0, 5, world, 5 ) <  0 );
    EXPECT_TRUE( string_view( hello ).compare( 6, 5, hello, 5 ) >  0 );
}

// 24.4.2.7 Searching:

#if nssv_HAVE_STARTS_WITH

TEST(StringViewTest, CanCheckForPrefixViewViaStartsWith)
{
    SCOPED_TRACE( "string_view: Allows to check for a prefix string_view via starts_with(), (1)" );

    char hello[] = "hello world";

    EXPECT_TRUE(     string_view( hello ).starts_with( string_view( hello ) ) );
    EXPECT_TRUE(     string_view( hello ).starts_with( string_view("hello") ) );
    EXPECT_FALSE( string_view( hello ).starts_with( string_view("world") ) );
}

TEST(StringViewTest, CanCheckForPrefixCharacterViaStartsWith)
{
    SCOPED_TRACE( "string_view: Allows to check for a prefix character via starts_with(), (2)" );

    char hello[] = "hello world";

    EXPECT_TRUE(     string_view( hello ).starts_with( 'h' ) );
    EXPECT_FALSE( string_view( hello ).starts_with( 'e' ) );
}

TEST(StringViewTest, CanCheckForPrefixCStringViaStartsWith)
{
    SCOPED_TRACE( "string_view: Allows to check for a prefix C-string via starts_with(), (3)" );

    char hello[] = "hello world";

    EXPECT_TRUE(     string_view( hello ).starts_with( hello ) );
    EXPECT_TRUE(     string_view( hello ).starts_with("hello") );
    EXPECT_FALSE( string_view( hello ).starts_with("world") );
}

#endif // nssv_HAVE_STARTS_WITH

#if nssv_HAVE_ENDS_WITH

TEST(StringViewTest, CanCheckForSuffixViewViaEndsWith)
{
    SCOPED_TRACE( "string_view: Allows to check for a suffix string_view via ends_with(), (1)" );

    char hello[] = "hello world";

    EXPECT_TRUE(     string_view( hello ).ends_with( string_view( hello ) ) );
    EXPECT_TRUE(     string_view( hello ).ends_with( string_view("world") ) );
    EXPECT_FALSE( string_view( hello ).ends_with( string_view("hello") ) );
}

TEST(StringViewTest, CanCheckForSuffixCharacterViaEndsWith)
{
    SCOPED_TRACE( "string_view: Allows to check for a suffix character via ends_with(), (2)" );

    char hello[] = "hello world";

    EXPECT_TRUE(     string_view( hello ).ends_with( 'd' ) );
    EXPECT_FALSE( string_view( hello ).ends_with( 'l' ) );
}

TEST(StringViewTest, CanCheckForSuffixCStringViaEndsWith)
{
    SCOPED_TRACE( "string_view: Allows to check for a suffix C-string via ends_with(), (3)" );

    char hello[] = "hello world";

    EXPECT_TRUE(     string_view( hello ).ends_with( hello ) );
    EXPECT_TRUE(     string_view( hello ).ends_with("world") );
    EXPECT_FALSE( string_view( hello ).ends_with("hello") );
}

#endif // nssv_HAVE_ENDS_WITH

TEST(StringViewTest, CanSearchForViewSubstrViaFind)
{
    SCOPED_TRACE( "string_view: Allows to search for a string_view substring, starting at position pos (default: 0) via find(), (1)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find( sv    ) == size_type( 0 ) );
    EXPECT_TRUE( sv.find( sv, 1 ) == string_view::npos );
    EXPECT_TRUE( sv.find( string_view("world" )    ) == size_type( 6 ) );
    EXPECT_TRUE( sv.find( string_view("world" ), 6 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.find( string_view("world" ), 7 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForCharacterViaFind)
{
    SCOPED_TRACE( "string_view: Allows to search for a character, starting at position pos (default: 0) via find(), (2)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find('h'    ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find('h', 1 ) == string_view::npos );
    EXPECT_TRUE( sv.find('w'    ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find('w', 6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find('w', 7 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForCStringSubstringViaFind)
{
    SCOPED_TRACE( "string_view: Allows to search for a C-string substring, starting at position pos and of length n via find(), (3)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find( hello , 0, sv.size() ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find( hello , 1, sv.size() ) == string_view::npos );
    EXPECT_TRUE( sv.find("world", 0, 5 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find("world", 6, 5 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find("world", 7, 4 ) == string_view::npos );
    EXPECT_TRUE( sv.find("world", 3, 0 ) == size_type( 3 )    );
}

TEST(StringViewTest, CanSearchForCStringSubstringViaFindWithDefaultPos)
{
    SCOPED_TRACE( "string_view: Allows to search for a C-string substring, starting at position pos (default: 0) via find(), (4)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find( hello     ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find( hello , 1 ) == string_view::npos );
    EXPECT_TRUE( sv.find("world"    ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find("world", 6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find("world", 7 ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForViewSubstrViaFind)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for a string_view substring, starting at position pos (default: npos) via rfind(), (1)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.rfind( sv    ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind( sv, 3 ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind( string_view(        )    ) == size_type(11 ) );
    EXPECT_TRUE( sv.rfind( string_view("world" )    ) == size_type( 6 ) );
    EXPECT_TRUE( sv.rfind( string_view("world" ), 6 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.rfind( string_view("world" ), 5 ) == string_view::npos );
    EXPECT_TRUE( sv.rfind( string_view("hello world, a longer text" ) ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForCharacterViaFind)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for a character, starting at position pos (default: npos) via rfind(), (2)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.rfind('h'    ) == size_type( 0 )    );
    EXPECT_TRUE( sv.rfind('e'    ) == size_type( 1 )    );
    EXPECT_TRUE( sv.rfind('e', 0 ) == string_view::npos );
    EXPECT_TRUE( sv.rfind('w'    ) == size_type( 6 )    );
    EXPECT_TRUE( sv.rfind('w', 6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.rfind('w', 5 ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForCStringSubstringViaFind)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for a C-string substring, starting at position pos and of length n via rfind(), (3)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.rfind( hello         ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind( hello ,  0, 5 ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind( hello ,  1, 5 ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind("world", 10, 5 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.rfind("world",  6, 5 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.rfind("world",  5, 5 ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForCStringSubstringViaFindWithDefaultPos)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for a C-string substring, starting at position pos (default: 0) via rfind(), (4)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.rfind( hello     ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind( hello , 3 ) == size_type( 0 ) );
    EXPECT_TRUE( sv.rfind("world"    ) == size_type( 6 ) );
    EXPECT_TRUE( sv.rfind("world", 6 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.rfind("world", 5 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForFirstOccurenceOfAnyCharacterInView)
{
    SCOPED_TRACE( "string_view: Allows to search for the first occurrence of any of the characters specified in a string view, starting at position pos (default: 0) via find_first_of(), (1)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_of( sv    ) == size_type( 0 ) );
    EXPECT_TRUE( sv.find_first_of( sv, 3 ) == size_type( 3 ) );
    EXPECT_TRUE( sv.find_first_of( string_view("xwo" )    ) == size_type( 4 ) );
    EXPECT_TRUE( sv.find_first_of( string_view("wdx" ), 6 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.find_first_of( string_view("wxy" ), 7 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForFirstOccurenceOfCharacter)
{
    SCOPED_TRACE( "string_view: Allows to search for a character, starting at position pos (default: 0) via find_first_of(), (2)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_of('h'    ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find_first_of('h', 1 ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_of('w'    ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_first_of('w', 6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_first_of('w', 7 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForFirstOccurenceOfCharactersInCStringInLenght)
{
    SCOPED_TRACE( "string_view: Allows to search for the first occurrence of any of the characters specified in a C-string, starting at position pos and of length n via find_first_of(), (3)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_of( hello , 0, sv.size() ) == size_type( 0 ) );
    EXPECT_TRUE( sv.find_first_of( hello , 1, sv.size() ) == size_type( 1 ) );
    EXPECT_TRUE( sv.find_first_of(  "xwy", 0, 3 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_first_of(  "xwy", 6, 3 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_first_of(  "xwy", 7, 3 ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_of(  "xyw", 0, 2 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForFirstOccurenceOfCharactersInCString)
{
    SCOPED_TRACE( "string_view: Allows to search for the first occurrence of any of the characters specified in a C-string, starting at position pos via find_first_of(), (4)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_of( hello , 0 ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find_first_of( hello , 1 ) == size_type( 1 )    );
    EXPECT_TRUE( sv.find_first_of(  "xwy", 0 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_first_of(  "xwy", 6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_first_of(  "xwy", 7 ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForLastOccurenceOfAnyCharacterInView)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for the last occurrence of any of the characters specified in a string view, starting at position pos (default: npos) via find_last_of(), (1)" );

    char hello[] = "hello world";
    char empty[] = "";
    string_view sv( hello );
    string_view sve( empty );

    EXPECT_TRUE( sv.find_last_of( sv    ) == size_type( 10 ) );
    EXPECT_TRUE( sv.find_last_of( sv, 3 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_of( string_view("xwo" )    ) == size_type( 7 ) );
    EXPECT_TRUE( sv.find_last_of( string_view("wdx" ), 6 ) == size_type( 6 ) );
    EXPECT_TRUE( sv.find_last_of( string_view("wxy" ), 7 ) == size_type( 6 ) );

    EXPECT_TRUE( sve.find_last_of( string_view("x") ) == string_view::npos );    // issue 20 (endless loop)
}

TEST(StringViewTest, CanBackwardsSearchForLastOccurenceOfCharacter)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for a character, starting at position pos (default: 0) via find_last_of(), (2)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_last_of('h'    ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find_last_of('l', 1 ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_of('w'    ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_last_of('w', 6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_last_of('w', 5 ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForLastOccurenceOfCharactersInCStringInLenght)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for the first occurrence of any of the characters specified in a C-string, starting at position pos and of length n via find_last_of(), (3)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_last_of( hello , 0, sv.size() ) == size_type( 0 ) );
    EXPECT_TRUE( sv.find_last_of( hello , 1, sv.size() ) == size_type( 1 ) );
    EXPECT_TRUE( sv.find_last_of("xwy", 10, 3 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_last_of("xwy",  6, 3 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_last_of("xwy",  5, 3 ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_of("xyw", 10, 2 ) == string_view::npos );
}

TEST(StringViewTest, CanBackwardsSearchForLastOccurenceOfCharactersInCString)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for the first occurrence of any of the characters specified in a C-string, starting at position pos via find_last_of(), (4)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_last_of( hello ,  0 ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find_last_of( hello ,  1 ) == size_type( 1 )    );
    EXPECT_TRUE( sv.find_last_of(  "xwy", 10 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_last_of(  "xwy",  6 ) == size_type( 6 )    );
    EXPECT_TRUE( sv.find_last_of(  "xwy",  5 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForFirstNotFoundCharacter)
{
    SCOPED_TRACE( "string_view: Allows to search for the first character not specified in a string view, starting at position pos (default: 0) via find_first_not_of(), (1)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_not_of( sv    ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_not_of( sv, 3 ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_not_of( string_view("helo "   )    ) == size_type(  6 ) );
    EXPECT_TRUE( sv.find_first_not_of( string_view("helo "   ), 6 ) == size_type(  6 ) );
    EXPECT_TRUE( sv.find_first_not_of( string_view("helo "   ), 7 ) == size_type(  8 ) );
    EXPECT_TRUE( sv.find_first_not_of( string_view("helo wr" )    ) == size_type( 10 ) );
}

TEST(StringViewTest, CanSearchForFirstNonMatchingCharacter)
{
    SCOPED_TRACE( "string_view: Allows to search for the first character not equal to the specified character, starting at position pos (default: 0) via find_first_not_of(), (2)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_not_of('h'     ) == size_type( 1 )    );
    EXPECT_TRUE( sv.find_first_not_of('h',  1 ) == size_type( 1 )    );
    EXPECT_TRUE( sv.find_first_not_of('w'     ) == size_type( 0 )    );
    EXPECT_TRUE( sv.find_first_not_of('w',  6 ) == size_type( 7 )    );
    EXPECT_TRUE( sv.find_first_not_of('d', 10 ) == string_view::npos );
}

TEST(StringViewTest, CanSearchForFirstNonEqualToAnyCharacterInCStringInLength)
{
    SCOPED_TRACE( "string_view: Allows to search for  the first character not equal to any of the characters specified in a C-string, starting at position pos and of length n via find_first_not_of(), (3)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_not_of( hello, 0, sv.size() ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_not_of( hello, 3, sv.size() ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_not_of( "helo "  , 0, 5     ) == size_type(  6 ) );
    EXPECT_TRUE( sv.find_first_not_of( "helo "  , 6, 5     ) == size_type(  6 ) );
    EXPECT_TRUE( sv.find_first_not_of( "helo "  , 7, 5     ) == size_type(  8 ) );
    EXPECT_TRUE( sv.find_first_not_of( "helo wr", 0, 7     ) == size_type( 10 ) );
    EXPECT_TRUE( sv.find_first_not_of( "he"     , 0, 1     ) == size_type(  1 ) );
}

TEST(StringViewTest, CanSearchForFirstNonEqualToAnyCharacterInCString)
{
    SCOPED_TRACE( "string_view: Allows to search for  the first character not equal to any of the characters specified in a C-string, starting at position pos via find_first_not_of(), (4)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_first_not_of( hello    , 0 ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_not_of( hello    , 3 ) == string_view::npos );
    EXPECT_TRUE( sv.find_first_not_of( "helo "  , 0 ) == size_type(  6 ) );
    EXPECT_TRUE( sv.find_first_not_of( "helo "  , 6 ) == size_type(  6 ) );
    EXPECT_TRUE( sv.find_first_not_of( "helo "  , 7 ) == size_type(  8 ) );
    EXPECT_TRUE( sv.find_first_not_of( "helo wr", 0 ) == size_type( 10 ) );
}

TEST(StringViewTest, CanBackwardsSearchForForstNonFoundCharacterInView)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for the first character not specified in a string view, starting at position pos (default: npos) via find_last_not_of(), (1)" );

    char hello[] = "hello world";
    char empty[] = "";
    string_view sv( hello );
    string_view sve( empty );

    EXPECT_TRUE( sv.find_last_not_of( sv    ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_not_of( sv, 3 ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_not_of( string_view("world " )    ) == size_type(  1 ) );
    EXPECT_TRUE( sv.find_last_not_of( string_view("heo "   ), 4 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_not_of( string_view("heo "   ), 3 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_not_of( string_view("heo "   ), 2 ) == size_type(  2 ) );
    EXPECT_TRUE( sv.find_last_not_of( string_view("x"      )    ) == size_type( 10 ) );

    EXPECT_TRUE( sve.find_last_not_of( string_view("x") ) == string_view::npos );    // issue 20 (endless loop)
}

TEST(StringViewTest, CanBackwardsSearchForFirstNonMatchingCharacter)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for the first character not equal to the specified character, starting at position pos (default: npos) via find_last_not_of(), (2)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_last_not_of('d'     ) == size_type( 9 ) );
    EXPECT_TRUE( sv.find_last_not_of('d', 10 ) == size_type( 9 ) );
    EXPECT_TRUE( sv.find_last_not_of('d',  9 ) == size_type( 9 ) );
    EXPECT_TRUE( sv.find_last_not_of('d',  8 ) == size_type( 8 ) );
    EXPECT_TRUE( sv.find_last_not_of('d',  0 ) == size_type( 0 ) );
}

TEST(StringViewTest, CanBackwardsSearchForFirstNonEqualToAnyCharacterInCStringInLength)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for  the first character not equal to any of the characters specified in a C-string, starting at position pos and of length n via find_last_not_of(), (3)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_last_not_of( hello, 0, sv.size() ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_not_of( hello, 3, sv.size() ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_not_of( "world ", 10, 6 ) == size_type(  1 ) );
    EXPECT_TRUE( sv.find_last_not_of( "heo "  ,  4, 4 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_not_of( "heo "  ,  3, 4 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_not_of( "heo "  ,  2, 4 ) == size_type(  2 ) );
    EXPECT_TRUE( sv.find_last_not_of( "x"             ) == size_type( 10 ) );
}

TEST(StringViewTest, CanBackwardsSearchForFirstNonEqualToAnyCharacterInCString)
{
    SCOPED_TRACE( "string_view: Allows to search backwards for  the first character not equal to any of the characters specified in a C-string, starting at position pos via find_last_not_of(), (4)" );

    char hello[] = "hello world";
    string_view sv( hello );

    EXPECT_TRUE( sv.find_last_not_of( hello    , 0 ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_not_of( hello    , 3 ) == string_view::npos );
    EXPECT_TRUE( sv.find_last_not_of( "world ", 10 ) == size_type(  1 ) );
    EXPECT_TRUE( sv.find_last_not_of( "heo "  ,  4 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_not_of( "heo "  ,  3 ) == size_type(  3 ) );
    EXPECT_TRUE( sv.find_last_not_of( "heo "  ,  2 ) == size_type(  2 ) );
    EXPECT_TRUE( sv.find_last_not_of( "x"          ) == size_type( 10 ) );
}

TEST(StringViewTest, CanCreateViewWithLiteralSV)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view, wstring_view, u16string_view, u32string_view via literal \"sv\"" );

#if nssv_CONFIG_STD_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_STD_DEFINED_LITERALS )
    using namespace nonstd::literals::string_view_literals;

    string_view sv1 =  "abc"sv;
    wstring_view sv2 = L"abc"sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 3 ) );

#if nssv_HAVE_WCHAR16_T
    u16string_view sv3 = u"abc"sv;
    EXPECT_TRUE( sv3.size() == size_type( 3 ) );
#endif
#if nssv_HAVE_WCHAR32_T
    u32string_view sv4 = U"abc"sv;
    EXPECT_TRUE( sv4.size() == size_type( 3 ) );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (nssv_CONFIG_STD_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteralSVInLiteralsStringViewLiteralsNamespace)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view via literal \"sv\", using namespace gmx::compat::literals::string_view_literals" );

#if nssv_CONFIG_STD_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_STD_DEFINED_LITERALS )
    using namespace gmx::compat::literals::string_view_literals;

    string_view sv1 = "abc\0\0def";
    string_view sv2 = "abc\0\0def"sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 8 ) );
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (nssv_CONFIG_STD_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteralSVInStringViewLiteralsNamespace)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view via literal \"sv\", using namespace gmx::compat::string_view_literals" );

#if nssv_CONFIG_STD_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_STD_DEFINED_LITERALS )
#if nssv_STD_SV_OR( nssv_HAVE_INLINE_NAMESPACE )
    using namespace gmx::compat::string_view_literals;

    string_view sv1 = "abc\0\0def";
    string_view sv2 = "abc\0\0def"sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 8 ) );
#else
    EXPECT_TRUE( !!"Inline namespaces for literal operator 'sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (nssv_CONFIG_STD_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteralSVInLiteralsNamespace)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view via literal \"sv\", using namespace gmx::compat::literals" );

#if nssv_CONFIG_STD_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_STD_DEFINED_LITERALS )
#if nssv_STD_SV_OR( nssv_HAVE_INLINE_NAMESPACE )
    using namespace gmx::compat::literals;

    string_view sv1 = "abc\0\0def";
    string_view sv2 = "abc\0\0def"sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 8 ) );
#else
    EXPECT_TRUE( !!"Inline namespaces for literal operator 'sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator 'sv' for string_view not available (nssv_CONFIG_STD_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteral_SV)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view, wstring_view, u16string_view, u32string_view via literal \"_sv\"" );

#if nssv_CONFIG_USR_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_USER_DEFINED_LITERALS )
    using namespace gmx::compat::literals::string_view_literals;

    string_view sv1 =  "abc"_sv;
    wstring_view sv2 = L"abc"_sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 3 ) );

#if nssv_HAVE_WCHAR16_T
    u16string_view sv3 = u"abc"_sv;
    EXPECT_TRUE( sv3.size() == size_type( 3 ) );
#endif
#if nssv_HAVE_WCHAR32_T
    u32string_view sv4 = U"abc"_sv;
    EXPECT_TRUE( sv4.size() == size_type( 3 ) );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (nssv_CONFIG_USR_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteral_SVInLiteralsStringViewLiteralsNamespace)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view via literal \"_sv\", using namespace gmx::compat::literals::string_view_literals" );

#if nssv_CONFIG_USR_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_USER_DEFINED_LITERALS )
    using namespace gmx::compat::literals::string_view_literals;

    string_view sv1 = "abc\0\0def";
    string_view sv2 = "abc\0\0def"_sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 8 ) );
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (nssv_CONFIG_USR_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteral_SVInStringViewLiteralsNamespace)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view via literal \"_sv\", using namespace gmx::compat::string_view_literals" );

#if nssv_CONFIG_USR_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_USER_DEFINED_LITERALS )
#if nssv_STD_SV_OR( nssv_HAVE_INLINE_NAMESPACE )
    using namespace gmx::compat::string_view_literals;

    string_view sv1 = "abc\0\0def";
    string_view sv2 = "abc\0\0def"_sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 8 ) );
#else
    EXPECT_TRUE( !!"Inline namespaces for literal operator '_sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (nssv_CONFIG_USR_SV_OPERATOR=0)." );
#endif
}

TEST(StringViewTest, CanCreateViewWithLiteral_SVInLiteralsNamespace)
{
    SCOPED_TRACE( "string_view: Allows to create a string_view via literal \"_sv\", using namespace gmx::compat::literals" );

#if nssv_CONFIG_USR_SV_OPERATOR
#if nssv_STD_SV_OR( nssv_HAVE_USER_DEFINED_LITERALS )
#if nssv_STD_SV_OR( nssv_HAVE_INLINE_NAMESPACE )
    using namespace gmx::compat::literals;

    string_view sv1 = "abc\0\0def";
    string_view sv2 = "abc\0\0def"_sv;

    EXPECT_TRUE( sv1.size() == size_type( 3 ) );
    EXPECT_TRUE( sv2.size() == size_type( 8 ) );
#else
    EXPECT_TRUE( !!"Inline namespaces for literal operator '_sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Literal operator '_sv' for string_view not available (nssv_CONFIG_USR_SV_OPERATOR=0)." );
#endif
}

// 24.4.3 Non-member comparison functions:

TEST(StringViewTest, CanCompareToViews)
{
    SCOPED_TRACE( "string_view: Allows to compare a string_view with another string_view" );

    char s[] = "hello";
    char t[] = "world";
    string_view sv( s );
    string_view tv( t );

    EXPECT_TRUE( sv.length() == size_type( 5 ) );
    EXPECT_TRUE( tv.length() == size_type( 5 ) );

    EXPECT_TRUE( sv == sv );
    EXPECT_TRUE( sv != tv );
    EXPECT_TRUE( sv <= sv );
    EXPECT_TRUE( sv <= tv );
    EXPECT_TRUE( sv <  tv );
    EXPECT_TRUE( tv >= tv );
    EXPECT_TRUE( tv >= sv );
    EXPECT_TRUE( tv >  sv );
}

TEST(StringViewTest, CanCompareViewToImplicitlyConvertedView)
{
    SCOPED_TRACE( "string_view: Allows to compare a string_view with an object with implicit conversion to string_view" );

#if nssv_CPP11_OR_GREATER
#if defined(_MSC_VER) && _MSC_VER != 1900
    char s[] = "hello";
    string_view sv( s );

    EXPECT_TRUE( sv == "hello"       );
    EXPECT_TRUE(       "hello" == sv );

    EXPECT_TRUE( sv != "world"       );
    EXPECT_TRUE(       "world" != sv );

    EXPECT_TRUE( sv <= "hello"       );
    EXPECT_TRUE(       "hello" <= sv );
    EXPECT_TRUE( sv <= "world"       );
    EXPECT_TRUE(       "aloha" <= sv );

    EXPECT_TRUE( sv <  "world"       );
    EXPECT_TRUE(       "aloha" <  sv );

    EXPECT_TRUE( sv >= "hello"       );
    EXPECT_TRUE(       "hello" >= sv );
    EXPECT_TRUE( sv >= "aloha"       );
    EXPECT_TRUE(       "world" >= sv );

    EXPECT_TRUE( sv >  "aloha"       );
    EXPECT_TRUE(       "world"  >  sv );
#else
    EXPECT_TRUE( !!"Comparison for types with implicit conversion to string_view not available (insufficient C++11 support of MSVC)." );
#endif
#else
    EXPECT_TRUE( !!"Comparison for types with implicit conversion to string_view not available (no C++11)." );
#endif
}

TEST(StringViewTest, EmptyViewsCompareAsEqual)
{
    SCOPED_TRACE( "string_view: Allows to compare empty string_view-s as equal" );

    string_view a, b;

    EXPECT_TRUE( a == b );
}

// 24.4.4 Inserters and extractors:

TEST(StringViewTest, CanPrintViewToPutputStream)
{
    SCOPED_TRACE ( "operator<<: Allows printing a string_view to an output stream" );

    std::ostringstream oss;
    char s[] = "hello";
    string_view sv( s );

    oss << sv << '\n'
        << std::right << std::setw(10) << sv << '\n'
        << sv << '\n'
        << std::setfill('.') << std::left << std::setw(10) << sv;

    EXPECT_TRUE( oss.str() == "hello\n     hello\nhello\nhello....." );
}

// 24.4.5 Hash support (C++11):

TEST(StringViewTest, HashOfViewIsEqualToHashOfString)
{
    SCOPED_TRACE ( "std::hash<>: Hash value of string_view equals hash value of corresponding string object" );

#if nssv_HAVE_STD_HASH
    EXPECT_TRUE( std::hash<string_view>()( "Hello, world!" ) == std::hash<std::string>()( "Hello, world!" ) );
#else
    EXPECT_TRUE( !!"std::hash is not available (no C++11)" );
#endif
}

TEST(StringViewTest, HashOfWStringViewIsEqualToHashOfString)
{
    SCOPED_TRACE ( "std::hash<>: Hash value of wstring_view equals hash value of corresponding string object" );

#if nssv_HAVE_STD_HASH
    EXPECT_TRUE( std::hash<wstring_view>()( L"Hello, world!" ) == std::hash<std::wstring>()( L"Hello, world!" ) );
#else
    EXPECT_TRUE( !!"std::hash is not available (no C++11)" );
#endif
}

TEST(StringViewTest, HashOfU16StringViewIsEqualToHashOfString)
{
    SCOPED_TRACE ( "std::hash<>: Hash value of u16string_view equals hash value of corresponding string object" );

#if nssv_HAVE_STD_HASH
#if nssv_HAVE_WCHAR16_T
#if nssv_HAVE_UNICODE_LITERALS
    EXPECT_TRUE( std::hash<u16string_view>()( u"Hello, world!" ) == std::hash<std::u16string>()( u"Hello, world!" ) );
#else
    EXPECT_TRUE( !!"Unicode literal u\"...\" is not available (no C++11)" );
#endif
#else
    EXPECT_TRUE( !!"std::u16string is not available (no C++11)" );
#endif
#else
    EXPECT_TRUE( !!"std::hash is not available (no C++11)" );
#endif
}

TEST(StringViewTest, HashOfU32StringViewIsEqualToHashOfString)
{
    SCOPED_TRACE ( "std::hash<>: Hash value of u32string_view equals hash value of corresponding string object" );

#if nssv_HAVE_STD_HASH
#if nssv_HAVE_WCHAR16_T
#if nssv_HAVE_UNICODE_LITERALS
    EXPECT_TRUE( std::hash<u32string_view>()( U"Hello, world!" ) == std::hash<std::u32string>()( U"Hello, world!" ) );
#else
    EXPECT_TRUE( !!"Unicode literal U\"...\" is not available (no C++11)" );
#endif
#else
    EXPECT_TRUE( !!"std::u32string is not available (no C++11)" );
#endif
#else
    EXPECT_TRUE( !!"std::hash is not available (no C++11)" );
#endif
}

// nonstd extension: conversions from and to std::basic_string

TEST(StringViewExtensionTest, CanConstructViewFromString)
{
    SCOPED_TRACE( "string_view: construct from std::string " "[extension]" );

#if nssv_USES_STD_STRING_VIEW
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_USES_STD_STRING_VIEW=1)." );
#elif nssv_CONFIG_CONVERSION_STD_STRING_CLASS_METHODS
    char hello[]  = "hello world";
    std::string s =  hello;

    string_view sv( hello );

    EXPECT_TRUE( sv.size() == s.size() );
    EXPECT_TRUE( sv.compare( s ) == 0  );
#else
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_CONFIG_CONVERSION_STD_STRING_CLASS_METHODS=0)." );
#endif
}

TEST(StringViewExtensionTest, CanConvertViewToStringViaExplicitOperator)
{
    SCOPED_TRACE( "string_view: convert to std::string via explicit operator " "[extension]" );

#if nssv_USES_STD_STRING_VIEW
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_USES_STD_STRING_VIEW=1)." );
#elif nssv_CONFIG_CONVERSION_STD_STRING_CLASS_METHODS
#if nssv_HAVE_EXPLICIT_CONVERSION
    char hello[] = "hello world";
    string_view sv( hello );

    std::string s( sv );
//  std::string t{ sv };

    EXPECT_TRUE( sv.size() == s.size() );
    EXPECT_TRUE( sv.compare( s ) == 0  );
#else
    EXPECT_TRUE( !!"explicit conversion is not available (no C++11)." );
#endif
#else
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_CONFIG_CONVERSION_STD_STRING_CLASS_METHODS=0)." );
#endif
}

TEST(StringViewExtensionTest, CanConvertViewToStringViaToString)
{
    SCOPED_TRACE( "string_view: convert to std::string via to_string() " "[extension]" );

#if nssv_USES_STD_STRING_VIEW
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_USES_STD_STRING_VIEW=1)." );
#elif nssv_CONFIG_CONVERSION_STD_STRING_CLASS_METHODS
    char hello[] = "hello world";
    string_view sv( hello );

    std::string s1 = sv.to_string();

    EXPECT_TRUE( sv.size() == s1.size() );
    EXPECT_TRUE( sv.compare( s1 ) == 0  );

    std::string s2 = sv.to_string( std::string::allocator_type() );

    EXPECT_TRUE( sv.size() == s2.size() );
    EXPECT_TRUE( sv.compare( s2 ) == 0  );
#else
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_CONFIG_CONVERSION_STD_STRING_CLASS_METHODS=0)." );
#endif
}

TEST(StringViewExtensionTest, CanConvertViewToStringViaToStringFreeFunction)
{
    SCOPED_TRACE( "to_string(): convert to std::string via to_string() " "[extension]" );

#if nssv_CONFIG_CONVERSION_STD_STRING_FREE_FUNCTIONS
    char hello[] = "hello world";
    string_view sv( hello );

    std::string s1 = to_string( sv );

    EXPECT_TRUE( sv.size() == s1.size() );
    EXPECT_TRUE( sv.compare( s1 ) == 0  );

    std::string s2 = to_string( sv, std::string::allocator_type() );

    EXPECT_TRUE( sv.size() == s2.size() );
    EXPECT_TRUE( sv.compare( s2 ) == 0  );

#else
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_CONFIG_CONVERSION_STD_STRING_FREE_FUNCTIONS=0)." );
#endif
}

TEST(StringViewExtensionTest, CanConvertViewToStringViewViaToStringView)
{
    SCOPED_TRACE( "to_string_view(): convert from std::string via to_string_view() " "[extension]" );

#if nssv_CONFIG_CONVERSION_STD_STRING_FREE_FUNCTIONS
    char hello[] = "hello world";
    std::string s = hello;

    string_view sv = to_string_view( s );

    EXPECT_TRUE( sv.size() == s.size() );
    EXPECT_TRUE( sv.compare( s ) == 0  );
#else
    EXPECT_TRUE( !!"Conversion to/from std::string is not available (nssv_CONFIG_CONVERSION_STD_STRING_FREE_FUNCTIONS=0)." );
#endif
}

} // anonymous namespace

} // namespace gmx

// GMX modification to suppress Doxygen checking
#endif // DOXYGEN
