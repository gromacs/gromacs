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
//
// Copyright 2014-2018 by Martin Moene
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// optional lite is inspired on std::optional by Fernando Cacciola and Andrzej Krzemienski
// and on expected lite by Martin Moene.

/*! \internal \file
 * \brief Tests for C++14-compatible implementation of std::optional.
 *
 * These tests are similer to those used in commit 7f4c6d4b4c8 of
 * https://github.com/martinmoene/optional-lite.git. The code has not
 * been linted with uncrustify so that any future updates to this
 * active repo can be incorporated easily, and //NOLINT comments added
 * to suppress clang-tidy warnings. The form of those
 * changes have been made to simplify the contents, while making it
 * easy to import any bug fixes that may appear in the source
 * repository.
 *
 * \todo Remove when requiring C++17, which has a standardized version
 * of std::optional.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_compat
 */
#include "gmxpre.h"

#include "gromacs/compat/optional.h"

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

#if optional_USES_STD_OPTIONAL && defined(__APPLE__)
# define opt_value( o ) *o
#else
# define opt_value( o )  o.value()
#endif

namespace {

struct nonpod { nonpod(){} };

struct Implicit { int x;          Implicit(int v) : x(v) {} };
struct Explicit { int x; explicit Explicit(int v) : x(v) {} };

bool operator==( Implicit a, Implicit b ) { return a.x == b.x; }
bool operator==( Explicit a, Explicit b ) { return a.x == b.x; }

std::ostream & operator<<( std::ostream & os, Implicit i ) { return os << "Implicit:" << i.x; }
std::ostream & operator<<( std::ostream & os, Explicit e ) { return os << "Explicit:" << e.x; }

// ensure comparison of pointers for lest:

// const void * lest_nullptr = 0;

// The following tracer code originates as Oracle from Optional by
// Andrzej Krzemienski, https://github.com/akrzemi1/Optional.

enum State
{
    /* 0 */ default_constructed,
    /* 1 */ value_copy_constructed,
    /* 2 */ value_move_constructed,
    /* 3 */ copy_constructed,
    /* 4 */ move_constructed,
    /* 5 */ move_assigned,
    /* 6 */ copy_assigned,
    /* 7 */ value_copy_assigned,
    /* 8 */ value_move_assigned,
    /* 9 */ moved_from,
    /*10 */ value_constructed
};

struct V
{
    State state;
    int   value;

    V(       ) : state( default_constructed ), value( deflt() ) {}
    V( int v ) : state( value_constructed   ), value( v       ) {}

    bool operator==( V const & rhs ) const { return state == rhs.state && value == rhs.value; }
    bool operator==( int       val ) const { return value == val; }

    static int deflt() { return 42; }
};

struct S
{
    State state;
    V     value;

    S(             ) : state( default_constructed    ) {}
    S( V const & v ) : state( value_copy_constructed ), value( v ) {}
    S( S const & s ) : state( copy_constructed       ), value( s.value        ) {}

    S & operator=( V const & v ) { state = value_copy_assigned; value = v; return *this; }
    S & operator=( const S & s ) { state = copy_assigned      ; value = s.value; return *this; }

#if optional_CPP11_OR_GREATER
    S(             V && v ) : state(  value_move_constructed ), value(  std::move( v       ) ) { v.state = moved_from; } //NOLINT(performance-move-const-arg) 
    S(             S && s ) : state(  move_constructed       ), value(  std::move( s.value ) ) { s.state = moved_from; } //NOLINT(performance-noexcept-move-constructor,performance-move-const-arg)

    S & operator=( V && v ) { state = value_move_assigned     ; value = std::move( v       ); v.state = moved_from; return *this; } //NOLINT(performance-noexcept-move-constructor,performance-move-const-arg,bugprone-use-after-move)
    S & operator=( S && s ) { state = move_assigned           ; value = std::move( s.value ); s.state = moved_from; return *this; } //NOLINT(performance-noexcept-move-constructor,performance-move-const-arg)
#endif

    bool operator==( S const & rhs ) const { return state == rhs.state && value == rhs.value; }
};

inline std::ostream & operator<<( std::ostream & os, V const & v )
{
    return os << "[V:" << toString( v.value ) << "]";
}

struct NoDefault
{
    NoDefault( NoDefault const & ) {} //NOLINT(readability-named-parameter)
    NoDefault & operator=( NoDefault const & ) { return *this; } //NOLINT(readability-named-parameter)

#if optional_CPP11_OR_GREATER
    NoDefault( NoDefault && ) = default;
    NoDefault & operator=( NoDefault && ) = default;
#endif

private:
    NoDefault();
};

struct CopyOnly
{
    CopyOnly( CopyOnly const & ) {} //NOLINT(readability-named-parameter)
    CopyOnly & operator=( CopyOnly const & ) { return *this; } //NOLINT(readability-named-parameter)

private:
    CopyOnly();
#if optional_CPP11_OR_GREATER
    CopyOnly( CopyOnly && ) = delete;
    CopyOnly & operator=( CopyOnly && ) = delete;
#endif
};

struct MoveOnly
{
#if optional_CPP11_OR_GREATER
    MoveOnly( MoveOnly && ) = default;
    MoveOnly & operator=( MoveOnly && ) = default;
#endif

private:
    MoveOnly();
    MoveOnly( MoveOnly const & );
    MoveOnly & operator=( MoveOnly const & );
};

struct NoDefaultCopyMove
{
    std::string text;
    NoDefaultCopyMove( std::string txt ) : text( txt ) {} //NOLINT(readability-named-parameter,performance-unnecessary-value-param)

private:
    NoDefaultCopyMove();
    NoDefaultCopyMove( NoDefaultCopyMove const & );
    NoDefaultCopyMove & operator=( NoDefaultCopyMove const & );
#if optional_CPP11_OR_GREATER
    NoDefaultCopyMove( NoDefaultCopyMove && ) = delete;
    NoDefaultCopyMove & operator=( NoDefaultCopyMove && ) = delete;
#endif
};

#if optional_CPP11_OR_GREATER
struct InitList
{
    std::vector<int> vec;
    char c;
    S s;

    InitList( std::initializer_list<int> il, char k, S const & t )
    : vec( il ), c( k ), s( t ) {}

    InitList( std::initializer_list<int> il, char k, S && t )
    : vec( il ), c( k ), s( std::move( t ) ) {}
};
#endif

} // anonymous namespace

//
// test specification:
//

TEST(OptionalTest, UnionCanContainNonPodTypes)
{
    SCOPED_TRACE("union: A C++03 union can only contain POD types");
    union U
    {
        char c;
#if optional_CPP11_OR_GREATER
        nonpod np;
#endif
    };
}

//
// optional member operations:
//

// construction:

TEST(OptionalTest, CanDefaultConstructEmpty)
{
    SCOPED_TRACE("optional: Allows to default construct an empty optional (1a)");
    optional<int> a;

    EXPECT_TRUE( !a );
}

TEST(OptionalTest, CanConstructFromNullopt)
{
    SCOPED_TRACE("optional: Allows to explicitly construct a disengaged, empty optional via nullopt (1b)");
    optional<int> a( nullopt );

    EXPECT_TRUE( !a );
}

TEST(OptionalTest, CanDefaultConstructUsingNonDefaultConstructibleType)
{
    SCOPED_TRACE("optional: Allows to default construct an empty optional with a non-default-constructible (1a)");
//  FAILS: NoDefault nd;
//  FAILS: NoDefaultCopyMove ndcm;
    optional<NoDefault> ond;
    optional<CopyOnly> oco;
    optional<MoveOnly> omo;
    optional<NoDefaultCopyMove> ondcm;

    EXPECT_TRUE( !ond );
    EXPECT_TRUE( !oco );
    EXPECT_TRUE( !omo );
    EXPECT_TRUE( !ondcm );
}

TEST(OptionalTest, CanCopyConstructFromEmptyOptional)
{
    SCOPED_TRACE("optional: Allows to copy-construct from empty optional (2)");
    optional<int> a;

    optional<int> b( a ); //NOLINT(performance-unnecessary-copy-initialization)

    EXPECT_TRUE( !b );
}

TEST(OptionalTest, CanMoveConstructFromEmptyOptional)
{
    SCOPED_TRACE("optional: Allows to move-construct from empty optional (C++11, 3)");
#if optional_CPP11_OR_GREATER
    optional<int> a;

    optional<int> b( std::move( a ) );

    EXPECT_TRUE( !b );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyConstructFromEmptyOptionalWithExplicitConversion)
{
    SCOPED_TRACE("optional: Allows to copy-construct from empty optional, explicit converting (C++11, 4a)");
#if optional_CPP11_OR_GREATER
    optional<int> a;

    optional<Explicit> b{ a };

    EXPECT_TRUE( !b );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyConstructFromEmptyOptionalNonExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to copy-construct from empty optional, non-explicit converting (4b)");
    optional<int> a;

    optional<Implicit> b( a );

    EXPECT_TRUE( !b );
}

TEST(OptionalTest, CanMoveConstructFromEmptyOptionalExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to move-construct from empty optional, explicit converting (C++11, 5a)");
#if optional_CPP11_OR_GREATER
    optional<int> a;

    optional<Explicit> b{ std::move( a ) };

    EXPECT_TRUE( !b );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveConstructFromEmptyOptionalNonExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to move-construct from empty optional, non-explicit converting (C++11, 5a)");
#if optional_CPP11_OR_GREATER
    optional<int> a;

    optional<Implicit> b( std::move( a ) );

    EXPECT_TRUE( !b );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyConstructFromNonEmptyOptional)
{
    SCOPED_TRACE("optional: Allows to copy-construct from non-empty optional (2)");
    optional<int> a( 7 );

    optional<int> b( a );

    EXPECT_TRUE(  b    );
    EXPECT_EQ  ( *b, 7 );
}

TEST(OptionalTest, CanCopyConstructFromNonEmptyOptionalExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to copy-construct from non-empty optional, explicit converting (C++11, 4a)");
#if optional_CPP11_OR_GREATER
    optional<int> a( 7 );

    optional<Explicit> b{ a };

    EXPECT_TRUE(  b                );
    EXPECT_EQ( *b, Explicit{7} );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyConstructFromNonEmptyOptionalNonExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to copy-construct from non-empty optional, non-explicit converting (4b)");
    optional<int> a( 7 );

    optional<Implicit> b( a );

    EXPECT_TRUE(  b                );
    EXPECT_EQ( *b, Implicit(7) );
}

TEST(OptionalTest, CanMoveConstructFromNonEmptyOptional)
{
    SCOPED_TRACE("optional: Allows to move-construct from non-empty optional (C++11, 3)");
#if optional_CPP11_OR_GREATER
    optional<int> a( 7 );

    optional<int> b( std::move( a ) );

    EXPECT_TRUE(  b      );
    EXPECT_EQ( *b, 7 );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveConstructFromNonEmptyOptionalExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to move-construct from non-empty optional, explicit converting (C++11, 5a)");
#if optional_CPP11_OR_GREATER
    optional<int> a( 7 );

    optional<Explicit> b{ std::move( a ) };

    EXPECT_TRUE(  b                );
    EXPECT_EQ( *b, Explicit{7} );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveConstructFromNonEmptyOptionalNonExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to move-construct from non-empty optional, non-explicit converting (C++11, 5b)");
#if optional_CPP11_OR_GREATER
    optional<int> a( 7 );

    optional<Implicit> b( std::move( a ) );

    EXPECT_TRUE(  b                );
    EXPECT_EQ( *b, Implicit(7) );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

namespace {

#if optional_CPP11_OR_GREATER
    void use_optional( nonstd::optional<Implicit> ) {} //NOLINT(readability-named-parameter,performance-unnecessary-value-param)
#else
    template< typename T >
    void use_optional( T ) {}
#endif

}

TEST(OptionalTest, CanCopyConstructFromLiteral)
{
    SCOPED_TRACE("optional: Allows to copy-construct from literal value (8)");
    use_optional( 7 );
    optional<int> a = 7;

    EXPECT_TRUE(  a      );
    EXPECT_EQ( *a, 7 );
}

TEST(OptionalTest, CanCopyConstructFromLiteralConverting)
{
    SCOPED_TRACE("optional: Allows to copy-construct from literal value, converting (8)");
    use_optional( '7' );
    optional<int> a = '7';

    EXPECT_TRUE(  a        );
    EXPECT_EQ( *a, '7' );
}

TEST(OptionalTest, CanCopyConstructFromValue)
{
    SCOPED_TRACE("optional: Allows to copy-construct from value (8)");
    const int i = 7;

    use_optional( i );
    optional<int> a( i );

    EXPECT_TRUE(  a      );
    EXPECT_EQ( *a, 7 );
}

TEST(OptionalTest, CanCopyConstructFromValueConverting)
{
    SCOPED_TRACE("optional: Allows to copy-construct from value, converting (8)");
    const char c = '7';

    use_optional( c );
    optional<int> a( c );

    EXPECT_TRUE(  a        );
    EXPECT_EQ( *a, '7' );
}

TEST(OptionalTest, CanMoveConstructFromValue)
{
    SCOPED_TRACE("optional: Allows to move-construct from value (C++11, 8b)");
#if optional_CPP11_OR_GREATER
    S s( 7 );

    optional<S> a( std::move( s ) );

    EXPECT_EQ( a->value, 7                );
    EXPECT_EQ( a->state, move_constructed );
    EXPECT_EQ(  s.state, moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveConstructFromValueExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to move-construct from value, explicit converting (C++11, 8a)");
#if optional_CPP11_OR_GREATER
    int seven = 7;

    optional<Explicit> a{ std::move( seven ) }; //NOLINT(performance-move-const-arg)

    EXPECT_TRUE(  a                );
    EXPECT_EQ( *a, Explicit{7} );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveConstructFromValueNonExplicitConverting)
{
    SCOPED_TRACE("optional: Allows to move-construct from value, non-explicit converting (C++11, 8b)");
#if optional_CPP11_OR_GREATER
    int seven = 7;
    optional<Implicit> a( std::move( seven ) ); //NOLINT(performance-move-const-arg)

    EXPECT_TRUE(  a                );
    EXPECT_EQ( *a, Implicit(7) );
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanInPlaceConstructFromLiteral)
{
    SCOPED_TRACE("optional: Allows to in-place construct from literal value (C++11, 6)");
#if optional_CPP11_OR_GREATER
    using pair_t = std::pair<char, int>;

    optional<pair_t> a( in_place, 'a', 7 );

    EXPECT_EQ( a->first , 'a' );
    EXPECT_EQ( a->second,  7  );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanInPlaceCopyConstructFromValue)
{
    SCOPED_TRACE("optional: Allows to in-place copy-construct from value (C++11, 6)");
#if optional_CPP11_OR_GREATER
    char c = 'a'; S s( 7 );
    using pair_t = std::pair<char, S>;

    optional<pair_t> a( in_place, c, s );

    EXPECT_EQ( a->first       , 'a' );
    EXPECT_EQ( a->second.value,  7  );
#if optional_USES_STD_OPTIONAL
    EXPECT_EQ( a->second.state, copy_constructed );
#else
    EXPECT_EQ( a->second.state, move_constructed );
#endif
    EXPECT_NE(         s.state , moved_from       );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanInPlaceMoveConstructFromValue)
{
    SCOPED_TRACE("optional: Allows to in-place move-construct from value (C++11, 6)");
#if optional_CPP11_OR_GREATER
    char c = 'a'; S s( 7 );
    using pair_t = std::pair<char, S>;

    optional<pair_t> a( in_place, c, std::move( s ) );

    EXPECT_EQ( a->first       , 'a' );
    EXPECT_EQ( a->second.value,  7  );
    EXPECT_EQ( a->second.state, move_constructed );
    EXPECT_EQ(         s.state, moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanInPlaceCopyConstructFromInitializerList)
{
    SCOPED_TRACE("optional: Allows to in-place copy-construct from initializer-list (C++11, 7)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    optional<InitList> a( in_place, { 7, 8, 9, }, 'a', s );

    EXPECT_EQ( a->vec[0] ,  7 );
    EXPECT_EQ( a->vec[1] ,  8 );
    EXPECT_EQ( a->vec[2] ,  9 );
    EXPECT_EQ( a->c      , 'a');
    EXPECT_EQ( a->s.value,  7 );
#if optional_USES_STD_OPTIONAL
    EXPECT_EQ( a->s.state , copy_constructed );
#else
    EXPECT_EQ( a->s.state , move_constructed );
#endif
    EXPECT_NE(    s.state , moved_from       );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanInPlaceMoveConstructFromInitializerList)
{
    SCOPED_TRACE("optional: Allows to in-place move-construct from initializer-list (C++11, 7)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    optional<InitList> a( in_place, { 7, 8, 9, }, 'a', std::move( s ) );

    EXPECT_EQ( a->vec[0]  ,  7  );
    EXPECT_EQ( a->vec[1]  ,  8  );
    EXPECT_EQ( a->vec[2]  ,  9  );
    EXPECT_EQ( a->c       , 'a' );
    EXPECT_EQ( a->s.value ,  7  );
    EXPECT_EQ( a->s.state , move_constructed );
    EXPECT_EQ(    s.state , moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

// assignment:

TEST(OptionalTest, CanAssignNulloptToDisengage)
{
    SCOPED_TRACE("optional: Allows to assign nullopt to disengage (1)");
    optional<int>  a( 7 );

    a = nullopt;

    EXPECT_TRUE( !a );
}

TEST(OptionalTest, CanCopyAssignBetweenEngagedAndDisengagedOptionals)
{
    SCOPED_TRACE("optional: Allows to copy-assign from/to engaged and disengaged optionals (2)");
        optional<int> d1;
        optional<int> d2;
        optional<int> e1( 123 );
        optional<int> e2( 987 );

    {
        SCOPED_TRACE("a disengaged optional assigned nullopt remains empty");
        d1 = nullopt;
        EXPECT_TRUE( !d1 );
    }
    {
        SCOPED_TRACE("a disengaged optional assigned an engaged optional obtains its value");
        d1 = e1;
        EXPECT_TRUE(  d1 );
        EXPECT_EQ( *d1 , 123 );
    }
    {
        SCOPED_TRACE("an engaged optional assigned an engaged optional obtains its value");
        e1 = e2;
        EXPECT_TRUE(  e1 );
        EXPECT_EQ( *e1 , 987 );
    }
    {
        SCOPED_TRACE("an engaged optional assigned nullopt becomes empty");
        e1 = nullopt;
        EXPECT_TRUE( !e1 );
    }
    {
        SCOPED_TRACE("a disengaged optional assigned a disengaged optional remains empty");
        d1 = d2;
        EXPECT_TRUE( !d1 );
    }
}

TEST(MakeOptionalTest, CanMoveAssignBetweenEngagedAndDisengagedOptionals)
{
    SCOPED_TRACE("optional: Allows to move-assign from/to engaged and disengaged optionals (C++11, 3)");
#if optional_CPP11_OR_GREATER
        optional<int> d1;
        optional<int> d2;
        optional<int> e1( 123 );
        optional<int> e2( 987 );

    {
        SCOPED_TRACE("a disengaged optional assigned nullopt remains empty");
        d1 = std::move( nullopt ); //NOLINT(performance-move-const-arg)
        EXPECT_TRUE( !d1 );
    }
    {
        SCOPED_TRACE("a disengaged optional assigned an engaged optional obtains its value");
        d1 = std::move( e1);
        EXPECT_TRUE(  d1 );
        EXPECT_EQ( *d1 , 123 );
    }
    {
        SCOPED_TRACE("an engaged optional assigned an engaged optional obtains its value");
        e1 = std::move( e2 );
        EXPECT_TRUE(  e1 );
        EXPECT_EQ( *e1 , 987 );
    }
    {
        SCOPED_TRACE("an engaged optional assigned nullopt becomes empty");
        e1 = std::move( nullopt ); //NOLINT(performance-move-const-arg)
        EXPECT_TRUE( !e1 );
    }
    {
        SCOPED_TRACE("a disengaged optional assigned a disengaged optional remains empty");
        d1 = std::move( d2);
        EXPECT_TRUE( !d1 );
    }
#else
    EXPECT_TRUE( !!"optional: move-assignment is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyAssignBetweenEngagedAndDisengagedOptionalsConverting)
{
    SCOPED_TRACE("optional: Allows to copy-assign from/to engaged and disengaged optionals, converting, (5)");
        optional<int>  d1;
        optional<char> d2;
        optional<int>  e1( 123 );
        optional<char> e2( '7' );

    {
        SCOPED_TRACE("a disengaged optional assigned an engaged optional obtains its value, converting");
        d1 = e2;
        EXPECT_TRUE(  d1 );
        EXPECT_EQ( *d1 , '7' );
    }
    {
        SCOPED_TRACE("an engaged optional assigned an engaged optional obtains its value, converting");
        e1 = e2;
        EXPECT_TRUE(  e1 );
        EXPECT_EQ( *e1 , '7' );
    }
    {
        SCOPED_TRACE("an engaged optional assigned a disengaged optional becomes empty, converting");
        e1 = d2;
        EXPECT_TRUE(  !e1 );
    }
    {
        SCOPED_TRACE("a disengaged optional assigned a disengaged optional remains empty, converting");
        d1 = d2;
        EXPECT_TRUE( !d1 );
    }
}

TEST(OptionalTest, CanMoveAssignBetweenEngagedAndDisengagedOptionalsConverting)
{
    SCOPED_TRACE("optional: Allows to move-assign from/to engaged and disengaged optionals, converting (C++11, 6)");
#if optional_CPP11_OR_GREATER
        optional<int>  d1;
        optional<char> d2;
        optional<int>  e1( 123 );
        optional<char> e2( '7' );

    {
        SCOPED_TRACE("a disengaged optional assigned an engaged optional obtains its value, converting");
        d1 = std::move( e2 );
        EXPECT_TRUE(  d1 );
        EXPECT_EQ( *d1 , '7' );
    }
    {
        SCOPED_TRACE("an engaged optional assigned an engaged optional obtains its value, converting");
        e1 = std::move( e2 ); //NOLINT(bugprone-use-after-move)
        EXPECT_TRUE(  e1 );
        EXPECT_EQ( *e1 , '7' );
    }
    {
        SCOPED_TRACE("an engaged optional assigned a disengaged optional becomes empty, converting");
        e1 = std::move( d2 );
        EXPECT_TRUE(  !e1 );
    }
    {
        SCOPED_TRACE("a disengaged optional assigned a disengaged optional remains empty, converting");
        d1 = std::move( d2 ); //NOLINT(bugprone-use-after-move)
        EXPECT_TRUE( !d1 );
    }
#else
    EXPECT_TRUE( !!"optional: move-assignment is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyAssignFromLiteral)
{
    SCOPED_TRACE("optional: Allows to copy-assign from literal value (4)");
    optional<int> a;

    a = 7;

    EXPECT_EQ( *a , 7 );
}

TEST(OptionalTest, CanCopyAssignFromValue)
{
    SCOPED_TRACE("optional: Allows to copy-assign from value (4)");
    const int i = 7;
    optional<int> a;

    a = i;

    EXPECT_EQ( *a , i );
}

TEST(OptionalTest, CanMoveAssignFromValue)
{
    SCOPED_TRACE("optional: Allows to move-assign from value (C++11, 4)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    optional<S> a;

    a = std::move( s );

    EXPECT_EQ( a->value , 7 );
    EXPECT_EQ( a->state , move_constructed );
    EXPECT_EQ(  s.state , moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: move-assign is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyEmplaceFromArguments)
{
    SCOPED_TRACE("optional: Allows to copy-emplace content from arguments (C++11, 7)");
#if optional_CPP11_OR_GREATER
    using pair_t = std::pair<char, S>;
    S s( 7 );
    optional<pair_t> a;

    a.emplace( 'a', s );

    EXPECT_EQ( a->first        , 'a' );
    EXPECT_EQ( a->second.value ,  7  );
    EXPECT_EQ( a->second.state , copy_constructed );
    EXPECT_NE(         s.state , moved_from       );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveEmplaceFromArguments)
{
    SCOPED_TRACE("optional: Allows to move-emplace content from arguments (C++11, 7)");
#if optional_CPP11_OR_GREATER
    using pair_t = std::pair<char, S>;
    S s( 7 );
    optional<pair_t> a;

    a.emplace( 'a', std::move( s ) );

    EXPECT_EQ( a->first        , 'a' );
    EXPECT_EQ( a->second.value ,  7  );
    EXPECT_EQ( a->second.state , move_constructed );
    EXPECT_EQ(         s.state , moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanCopyEmplaceFromInitializerListAndArguments)
{
    SCOPED_TRACE("optional: Allows to copy-emplace content from intializer-list and arguments (C++11, 8)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    optional<InitList> a;

    a.emplace( { 7, 8, 9, }, 'a', s );

    EXPECT_EQ( a->vec[0]  ,  7  );
    EXPECT_EQ( a->vec[1]  ,  8  );
    EXPECT_EQ( a->vec[2]  ,  9  );
    EXPECT_EQ( a->c       , 'a' );
    EXPECT_EQ( a->s.value ,  7  );
    EXPECT_EQ( a->s.state , copy_constructed );
    EXPECT_NE(    s.state , moved_from       );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanMoveEmplaceFromInitializerListAndArguments)
{
    SCOPED_TRACE("optional: Allows to move-emplace content from intializer-list and arguments (C++11, 8)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    optional<InitList> a;

    a.emplace( { 7, 8, 9, }, 'a', std::move( s ) );

    EXPECT_EQ( a->vec[0]  ,  7  );
    EXPECT_EQ( a->vec[1]  ,  8  );
    EXPECT_EQ( a->vec[2]  ,  9  );
    EXPECT_EQ( a->c       , 'a' );
    EXPECT_EQ( a->s.value ,  7               );
    EXPECT_EQ( a->s.state , move_constructed );
    EXPECT_EQ(    s.state , moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

// swap:

class OptionalMemberSwapTest : public ::testing::Test
{
    public:
        optional<int> d1;
        optional<int> d2;
        optional<int> e1{ 42 };
        optional<int> e2{ 7 };
};

TEST_F(OptionalMemberSwapTest, CanSwapDisengagedWithDisengaged)
    {
        SCOPED_TRACE("swap disengaged with disengaged optional");
        d1.swap( d2 );
        EXPECT_TRUE( !d1 );
    }

TEST_F(OptionalMemberSwapTest, CanSwapEngagedWithDisengaged)
    {
        SCOPED_TRACE("swap engaged with engaged optional");
        e1.swap( e2 );
        EXPECT_TRUE(  e1  );
        EXPECT_TRUE(  e2 );
        EXPECT_EQ( *e1 , 7 );
        EXPECT_EQ( *e2 , 42 );
    }
TEST_F(OptionalMemberSwapTest, CanSwapDisengagedWithEngaged)
    {
        SCOPED_TRACE("swap disengaged with engaged optional");
        d1.swap( e1 );
        EXPECT_TRUE(  d1 );
        EXPECT_TRUE( !e1 );
        EXPECT_EQ( *d1 , 42 );
    }
TEST_F(OptionalMemberSwapTest, CanSwapEngagedWithEngaged)
    {
        SCOPED_TRACE("swap engaged with disengaged optional");
        e1.swap( d1 );
        EXPECT_TRUE(  d1 );
        EXPECT_TRUE( !e1 );
        EXPECT_EQ( *d1 , 42 );
    }

// observers:

class OptionalImplicitValueTest : public ::testing::Test
{
    public:
        optional<Implicit>        e{ Implicit( 42 ) };
        optional<Implicit> const ce{ Implicit( 42 ) };
};

TEST_F(OptionalImplicitValueTest, CanObtainValueConst)
    {
        SCOPED_TRACE("operator->() yields pointer to value (const)");
        EXPECT_EQ(  ce->x , 42 );
    }
TEST_F(OptionalImplicitValueTest, CanObtainValueNonConst)
    {
        SCOPED_TRACE("operator->() yields pointer to value (non-const)");
        e->x = 7;
        EXPECT_EQ(  e->x , 7 );
    }

TEST_F(OptionalImplicitValueTest, CanObtainMovedFromValueConst)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("operator->() yields pointer to value (const)");
        EXPECT_EQ(  std::move( ce )->x , 42 ); //NOLINT(performance-move-const-arg)
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }
TEST_F(OptionalImplicitValueTest, CanObtainMovedFromValueNonConst)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("operator->() yields pointer to value (non-const)");
        e->x = 7;
        EXPECT_EQ(  std::move( e )->x , 7 );
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }

class OptionalIntValueTest : public ::testing::Test
{
    public:
        optional<int> a;
        optional<int> b{ 7 };
        optional<int> d;
        optional<int> const cd;
        optional<int>        e{ 42 };
        optional<int> const ce{ 42 };
        optional<int> d1;
        optional<int> d2;
        optional<int> e1{ 42 };
        optional<int> e2{ 7 };
};

TEST_F(OptionalIntValueTest, CanObtainValueFromDereferenceOperatorConst)
    {
        SCOPED_TRACE("operator*() yields value (const)");
        EXPECT_EQ( *ce , 42 );
    }
TEST_F(OptionalIntValueTest, CanObtainValueFromDereferenceOperatorNonConst)
    {
        SCOPED_TRACE("operator*() yields value (non-const)");
        *e = 7;
        EXPECT_EQ( *e , 7 );
    }

TEST_F(OptionalIntValueTest, CanObtainMovedValueFromDereferenceOperatorConst)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("operator*() yields value (const)");
        EXPECT_EQ( *(std::move( ce )) , 42 ); //NOLINT(performance-move-const-arg)
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }
TEST_F(OptionalIntValueTest, CanObtainMovedValueFromDereferenceOperatorNonConst)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("operator*() yields value (non-const)");
        *e = 7;
        EXPECT_EQ( *(std::move( e )) , 7 );
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }

TEST_F(OptionalIntValueTest, CanObtainHasValueViaOperatorBool)
{
    SCOPED_TRACE("optional: Allows to obtain has_value() via operator bool()");

    EXPECT_FALSE( a );
    EXPECT_TRUE(     b );
}

TEST_F(OptionalIntValueTest, CanObtainValueViaValueMethodConst)
    {
        SCOPED_TRACE("value() yields value (const)");
        EXPECT_EQ( opt_value( ce ) , 42 );
    }
TEST_F(OptionalIntValueTest, CanObtainValueViaValueMethodNonConst)
    {
        SCOPED_TRACE("value() yields value (non-const)");
        EXPECT_EQ( opt_value( e ) , 42 );
    }

TEST_F(OptionalIntValueTest, CanObtainMovedValueViaValueMethodConst)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("value() yields value (const)");
        EXPECT_EQ( opt_value( std::move( ce ) ) , 42 ); //NOLINT(performance-move-const-arg)
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }
TEST_F(OptionalIntValueTest, CanObtainMovedValueViaValueMethodNonConst)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("value() yields value (non-const)");
        EXPECT_EQ( opt_value( std::move( e ) ) , 42 );
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }

TEST_F(OptionalIntValueTest, CanObtainValueFromNonEmptyOptionalViaValueOrMethod)
    {
        SCOPED_TRACE("value_or( 7 ) yields value for non-empty optional");
        EXPECT_EQ( e.value_or( 7 ) , 42 );
    }
TEST_F(OptionalIntValueTest, CanObtainDefaultFromEmptyOptionalViaValueOrMethod)
    {
        SCOPED_TRACE("value_or( 7 ) yields default for empty optional");
        EXPECT_EQ( d.value_or( 7 ) , 7 );
    }

TEST_F(OptionalIntValueTest, CanObtainMovedFromValueForLValuesViaValueOrMethod)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("for l-values");
        EXPECT_EQ( d.value_or( 7 ) ,  7 );
        EXPECT_EQ( e.value_or( 7 ) , 42 );
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }
TEST_F(OptionalIntValueTest, CanObtainMovedFromValueForRValuesViaValueOrMethod)
    {
#if optional_CPP11_OR_GREATER
        SCOPED_TRACE("for r-values");
        EXPECT_EQ( std::move( d ).value_or( 7 ) ,  7 );
        EXPECT_EQ( std::move( e ).value_or( 7 ) , 42 );
#else
        EXPECT_TRUE( !!"optional: move-semantics are not available (no C++11)" );
#endif
    }

TEST_F(OptionalIntValueTest, ThrowsBadOptionalAccessAtDisengagedAccessForLValues)
{
    SCOPED_TRACE("optional: Throws bad_optional_access at disengaged access");
#if optional_USES_STD_OPTIONAL && defined(__APPLE__)
    EXPECT_TRUE( true );
#else

    {
        SCOPED_TRACE("for l-values");
        EXPECT_THROW(  d.value(), bad_optional_access );
        EXPECT_THROW( cd.value(), bad_optional_access );
    }
#endif
}
TEST_F(OptionalIntValueTest, ThrowsBadOptionalAccessAtDisengagedAccessForRValues)
{
    SCOPED_TRACE("optional: Throws bad_optional_access at disengaged access");
#if optional_USES_STD_OPTIONAL && defined(__APPLE__)
    EXPECT_TRUE( true );
#else

# if optional_CPP11_OR_GREATER
    {
        SCOPED_TRACE("for r-values");
        EXPECT_THROW( std::move(  d ).value(), bad_optional_access );
        EXPECT_THROW( std::move( cd ).value(), bad_optional_access ); //NOLINT(performance-move-const-arg)
    }
# endif
#endif
}

TEST_F(OptionalIntValueTest, ThrowsBadOptionalAccessWithNonEmptyWhatMethod)
{
    SCOPED_TRACE("optional: Throws bad_optional_access with non-empty what()");
    try
    {
        (void) d.value();
    }
    catch( bad_optional_access const & e )
    {
        EXPECT_TRUE( ! std::string( e.what() ).empty() );
    }
}

// modifiers:

TEST(OptionalTest, CanResetContent)
{
    SCOPED_TRACE("optional: Allows to reset content");
    optional<int> a = 7;

    a.reset();

    EXPECT_FALSE( a.has_value() );
}

//
// optional non-member functions:
//

TEST_F(OptionalIntValueTest, CanNonMemberSwapDisengagnedWithDisengaged)
    {
        SCOPED_TRACE("swap disengaged with disengaged optional");
        swap( d1, d2 );
        EXPECT_TRUE( !d1 );
    }
TEST_F(OptionalIntValueTest, CanNonMemberSwapEngagnedWithEngaged)
    {
        SCOPED_TRACE("swap engaged with engaged optional");
        swap( e1, e2 );
        EXPECT_TRUE(  e1  );
        EXPECT_TRUE(  e2 );
        EXPECT_EQ( *e1 , 7 );
        EXPECT_EQ( *e2 , 42 );
    }
TEST_F(OptionalIntValueTest, CanNonMemberSwapDisengagnedWithEngaged)
    {
        SCOPED_TRACE("swap disengaged with engaged optional");
        swap( d1, e1 );
        EXPECT_TRUE(  d1 );
        EXPECT_TRUE( !e1 );
        EXPECT_EQ( *d1 , 42 );
    }
TEST_F(OptionalIntValueTest, CanNonMemberSwapEngagnedWithDisengaged)
    {
        SCOPED_TRACE("swap engaged with disengaged optional");
        swap( e1, d1 );
        EXPECT_TRUE(  d1 );
        EXPECT_TRUE( !e1 );
        EXPECT_EQ( *d1 , 42 );
    }

template< typename R, typename S, typename T >
void relop()
{
        optional<R> d;
        optional<S> e1( 6 );
        optional<T> e2( 7 );

    { SCOPED_TRACE( "engaged    == engaged"    ); EXPECT_TRUE(   e1 == e1  ); }
    { SCOPED_TRACE( "engaged    == disengaged" ); EXPECT_TRUE( !(e1 == d ) ); }
    { SCOPED_TRACE( "disengaged == engaged"    ); EXPECT_TRUE( !(d  == e1) ); }

    { SCOPED_TRACE( "engaged    != engaged"    ); EXPECT_TRUE(   e1 != e2  ); }
    { SCOPED_TRACE( "engaged    != disengaged" ); EXPECT_TRUE(   e1 != d   ); }
    { SCOPED_TRACE( "disengaged != engaged"    ); EXPECT_TRUE(   d  != e2  ); }

    { SCOPED_TRACE( "engaged    <  engaged"    ); EXPECT_TRUE(   e1 <  e2  ); }
    { SCOPED_TRACE( "engaged    <  disengaged" ); EXPECT_TRUE( !(e1 <  d ) ); }
    { SCOPED_TRACE( "disengaged <  engaged"    ); EXPECT_TRUE(   d  <  e2  ); }

    { SCOPED_TRACE( "engaged    <= engaged"    ); EXPECT_TRUE(   e1 <= e1  ); }
    { SCOPED_TRACE( "engaged    <= engaged"    ); EXPECT_TRUE(   e1 <= e2  ); }
    { SCOPED_TRACE( "engaged    <= disengaged" ); EXPECT_TRUE( !(e1 <= d ) ); }
    { SCOPED_TRACE( "disengaged <= engaged"    ); EXPECT_TRUE(   d  <= e2  ); }

    { SCOPED_TRACE( "engaged    >  engaged"    ); EXPECT_TRUE(   e2 >  e1  ); }
    { SCOPED_TRACE( "engaged    >  disengaged" ); EXPECT_TRUE(   e2 >  d   ); }
    { SCOPED_TRACE( "disengaged >  engaged"    ); EXPECT_TRUE( !(d  >  e1) ); }

    { SCOPED_TRACE( "engaged    >= engaged"    ); EXPECT_TRUE(   e1 >= e1  ); }
    { SCOPED_TRACE( "engaged    >= engaged"    ); EXPECT_TRUE(   e2 >= e1  ); }
    { SCOPED_TRACE( "engaged    >= disengaged" ); EXPECT_TRUE(   e2 >= d   ); }
    { SCOPED_TRACE( "disengaged >= engaged"    ); EXPECT_TRUE( !(d  >= e1) ); }

    { SCOPED_TRACE( "disengaged == nullopt"    ); EXPECT_TRUE(  (d       == nullopt) ); }
    { SCOPED_TRACE( "nullopt    == disengaged" ); EXPECT_TRUE(  (nullopt == d      ) ); }
    { SCOPED_TRACE( "engaged    == nullopt"    ); EXPECT_TRUE(  (e1      != nullopt) ); }
    { SCOPED_TRACE( "nullopt    == engaged"    ); EXPECT_TRUE(  (nullopt != e1     ) ); }
    { SCOPED_TRACE( "disengaged == nullopt"    ); EXPECT_TRUE( !(d       <  nullopt) ); }
    { SCOPED_TRACE( "nullopt    == disengaged" ); EXPECT_TRUE( !(nullopt <  d      ) ); }
    { SCOPED_TRACE( "disengaged == nullopt"    ); EXPECT_TRUE(  (d       <= nullopt) ); }
    { SCOPED_TRACE( "nullopt    == disengaged" ); EXPECT_TRUE(  (nullopt <= d      ) ); }
    { SCOPED_TRACE( "disengaged == nullopt"    ); EXPECT_TRUE( !(d       >  nullopt) ); }
    { SCOPED_TRACE( "nullopt    == disengaged" ); EXPECT_TRUE( !(nullopt >  d      ) ); }
    { SCOPED_TRACE( "disengaged == nullopt"    ); EXPECT_TRUE(  (d       >= nullopt) ); }
    { SCOPED_TRACE( "nullopt    == disengaged" ); EXPECT_TRUE(  (nullopt >= d      ) ); }

    { SCOPED_TRACE( "engaged    == value"      ); EXPECT_TRUE( e1 == 6  ); }
    { SCOPED_TRACE( "value     == engaged"     ); EXPECT_TRUE(  6 == e1 ); }
    { SCOPED_TRACE( "engaged   != value"       ); EXPECT_TRUE( e1 != 7  ); }
    { SCOPED_TRACE( "value     != engaged"     ); EXPECT_TRUE(  6 != e2 ); }
    { SCOPED_TRACE( "engaged   <  value"       ); EXPECT_TRUE( e1 <  7  ); }
    { SCOPED_TRACE( "value     <  engaged"     ); EXPECT_TRUE(  6 <  e2 ); }
    { SCOPED_TRACE( "engaged   <= value"       ); EXPECT_TRUE( e1 <= 7  ); }
    { SCOPED_TRACE( "value     <= engaged"     ); EXPECT_TRUE(  6 <= e2 ); }
    { SCOPED_TRACE( "engaged   >  value"       ); EXPECT_TRUE( e2 >  6  ); }
    { SCOPED_TRACE( "value     >  engaged"     ); EXPECT_TRUE(  7 >  e1 ); }
    { SCOPED_TRACE( "engaged   >= value"       ); EXPECT_TRUE( e2 >= 6  ); }
    { SCOPED_TRACE( "value     >= engaged"     ); EXPECT_TRUE(  7 >= e1 ); }
}

TEST(OptionalTest, ProvidesRelationalOperators)
{
    SCOPED_TRACE("optional: Provides relational operators");
    relop<int, int, int>();
}

TEST(OptionalTest, ProvidesMixedTypeRelationalOperators)
{
    SCOPED_TRACE("optional: Provides mixed-type relational operators");
    relop<char, int, long>();
}

TEST(MakeOptionalTest, CanCopyConstruct)
{
    SCOPED_TRACE("make_optional: Allows to copy-construct optional");
    S s( 7 );

    EXPECT_EQ( make_optional( s )->value , 7          );
    EXPECT_NE(                   s.state , moved_from );
}

TEST(MakeOptionalTest, CanMoveConstruct)
{
    SCOPED_TRACE("make_optional: Allows to move-construct optional (C++11)");
#if optional_CPP11_OR_GREATER
    S s( 7 );

    EXPECT_EQ( make_optional( std::move( s ) )->value , 7          );
    EXPECT_EQ(                                s.state , moved_from ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: move-construction is not available (no C++11)" );
#endif
}

TEST(MakeOptionalTest, CanInPlaceCopyConstructFromArguments)
{
    SCOPED_TRACE("make_optional: Allows to in-place copy-construct optional from arguments (C++11)");
#if optional_CPP11_OR_GREATER
    using pair_t = std::pair<char, S>;

    S s( 7);
    auto a = make_optional<pair_t>( 'a', s );

    EXPECT_EQ( a->first        , 'a' );
    EXPECT_EQ( a->second.value ,  7  );
#if optional_USES_STD_OPTIONAL
    EXPECT_EQ( a->second.state , copy_constructed );
#else
    EXPECT_EQ( a->second.state , move_constructed );
#endif
    EXPECT_NE(         s.state , moved_from       );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(MakeOptionalTest, CanInPlaceMoveConstructFromArguments)
{
    SCOPED_TRACE("make_optional: Allows to in-place move-construct optional from arguments (C++11)");
#if optional_CPP11_OR_GREATER
    using pair_t = std::pair<char, S>;

    S s( 7 );
    auto a = make_optional<pair_t>( 'a', std::move( s ) );

    EXPECT_EQ( a->first        , 'a' );
    EXPECT_EQ( a->second.value ,  7  );
    EXPECT_EQ( a->second.state , move_constructed );
    EXPECT_EQ(         s.state , moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(MakeOptionalTest, CanInPlaceCopyConstructFromInitializerListAndArguments)
{
    SCOPED_TRACE("make_optional: Allows to in-place copy-construct optional from initializer-list and arguments (C++11)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    auto a = make_optional<InitList>( { 7, 8, 9, }, 'a', s );

    EXPECT_EQ( a->vec[0]  ,  7  );
    EXPECT_EQ( a->vec[1]  ,  8  );
    EXPECT_EQ( a->vec[2]  ,  9  );
    EXPECT_EQ( a->c       , 'a' );
    EXPECT_EQ( a->s.value ,  7  );
#if optional_USES_STD_OPTIONAL
    EXPECT_EQ( a->s.state , copy_constructed );
#else
    EXPECT_EQ( a->s.state , move_constructed );
#endif
    EXPECT_NE(    s.state , moved_from       );
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(MakeOptionalTest, CanInPlaceMoveConstructFromInitializerListAndArguments)
{
    SCOPED_TRACE("make_optional: Allows to in-place move-construct optional from initializer-list and arguments (C++11)");
#if optional_CPP11_OR_GREATER
    S s( 7 );
    auto a = make_optional<InitList>( { 7, 8, 9, }, 'a', std::move( s ) );

    EXPECT_EQ( a->vec[0]  ,  7  );
    EXPECT_EQ( a->vec[1]  ,  8  );
    EXPECT_EQ( a->vec[2]  ,  9  );
    EXPECT_EQ( a->c       , 'a' );
    EXPECT_EQ( a->s.value ,  7  );
    EXPECT_EQ( a->s.state , move_constructed );
    EXPECT_EQ(    s.state , moved_from       ); //NOLINT(bugprone-use-after-move)
#else
    EXPECT_TRUE( !!"optional: in-place construction is not available (no C++11)" );
#endif
}

TEST(OptionalTest, CanProduceHash)
{
    SCOPED_TRACE("std::hash<>: Allows to obtain hash (C++11)");
#if optional_CPP11_OR_GREATER
    const auto a = optional<int>( 7 );
    const auto b = optional<int>( 7 );

    EXPECT_EQ( std::hash<optional<int>>{}( a ) , std::hash<optional<int>>{}( b ) );
#else
    EXPECT_GT( !!"std::hash<>: std::hash<, is not available (no C++11)" );
#endif
}


//
// Negative tests:
//

//
// Tests that print information:
//

struct Struct{ Struct(){} };

#if !defined(optional_FEATURE_MAX_ALIGN_HACK) || !optional_FEATURE_MAX_ALIGN_HACK

#define optional_OUTPUT_ALIGNMENT_OF( type ) \
    "alignment_of<" #type ">: " <<  \
     alignment_of<type>::value  << "\n" <<

TEST(OptionalTest, ShowAlignmentDependingOnBaseType)
{
    SCOPED_TRACE("alignment_of: Show alignment of various types");

#if optional_CPP11_OR_GREATER
    using std::alignment_of;
#else
    using ::nonstd::optional_lite::detail::alignment_of;
#endif

    std::cout <<
        optional_OUTPUT_ALIGNMENT_OF( char )
        optional_OUTPUT_ALIGNMENT_OF( short )
        optional_OUTPUT_ALIGNMENT_OF( int )
        optional_OUTPUT_ALIGNMENT_OF( long )
        optional_OUTPUT_ALIGNMENT_OF( float )
        optional_OUTPUT_ALIGNMENT_OF( double )
        optional_OUTPUT_ALIGNMENT_OF( long double )
        optional_OUTPUT_ALIGNMENT_OF( Struct )
         "";
}
#undef optional_OUTPUT_ALIGNMENT_OF
#endif

#define optional_OUTPUT_SIZEOF( type ) \
    "sizeof( optional<" #type "> ): " << \
     sizeof( optional<   type>   )    << " (" << sizeof(type) << ")\n" <<

TEST(OptionalTest, ShowSizeDependingOnBaseType)
{
    std::cout <<
#if !optional_USES_STD_OPTIONAL
        "sizeof( nonstd::optional_lite::detail::storage_t<char> ): " <<
         sizeof( nonstd::optional_lite::detail::storage_t<char> )    << "\n" <<
#endif
         optional_OUTPUT_SIZEOF( char )
         optional_OUTPUT_SIZEOF( short )
         optional_OUTPUT_SIZEOF( int )
         optional_OUTPUT_SIZEOF( long )
         optional_OUTPUT_SIZEOF( float )
         optional_OUTPUT_SIZEOF( double )
         optional_OUTPUT_SIZEOF( long double )
         optional_OUTPUT_SIZEOF( Struct )
         "";
}
#undef optional_OUTPUT_SIZEOF

// GMX modification to suppress Doxygen checking
#endif // DOXYGEN

} // namespace gmx
