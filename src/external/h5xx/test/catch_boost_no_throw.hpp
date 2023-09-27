/*
 * Copyright © 2013  Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

//
// The following code snippet has been taken from
// http://stackoverflow.com/questions/15133259/boost-check-no-throw-how-to-get-exception-message-printed
//
// 2016/02: In Boost 1.60 the macro BOOST_CHECK_IMPL was removed, so the
// redefine below does not work any more and was put inside a version check.
// TODO : implement version independent solution
//


#ifndef CATCH_BOOST_NO_THROW_HPP
#define CATCH_BOOST_NO_THROW_HPP

#include <boost/test/unit_test.hpp>

#include <boost/preprocessor/stringize.hpp>
#include <exception>
#include <iostream>

// --- macro redefine only works with Boost versions before 1.60
#if BOOST_VERSION < 106000

#ifdef BOOST_CHECK_NO_THROW_IMPL
#   undef BOOST_CHECK_NO_THROW_IMPL
#endif
#define BOOST_CHECK_NO_THROW_IMPL( S, TL )                                                      \
    try {                                                                                       \
    S;                                                                                          \
    BOOST_CHECK_IMPL( true, "no exceptions thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG ); } \
    catch( const std::exception & e ) {                                                         \
    std::cerr << std::endl                                                                      \
    << "-----------------------------------------------" << std::endl                   \
    << std::endl << "exception message: " << e.what() << std::endl;                 \
    BOOST_CHECK_IMPL( false, "exception thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG );      \
    }                                                                                           \
    catch( ... ) {                                                                              \
    std::cerr << std::endl                                                                      \
    << "-----------------------------------------------" << std::endl                   \
    << std::endl << "exception message : <unkown exception>" << std::endl;          \
    BOOST_CHECK_IMPL( false, "exception thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG );      \
    }                                                                                           \
    /**/

#endif  // BOOST_VERSION

#define BOOST_WARN_NO_THROW( S )            BOOST_CHECK_NO_THROW_IMPL( S, WARN )
#define BOOST_CHECK_NO_THROW( S )           BOOST_CHECK_NO_THROW_IMPL( S, CHECK )
#define BOOST_REQUIRE_NO_THROW( S )         BOOST_CHECK_NO_THROW_IMPL( S, REQUIRE )

#endif // ! CATCH_BOOST_NO_THROW_HPP
