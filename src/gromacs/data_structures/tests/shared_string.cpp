/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015 by the GROMACS development team, led by
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
/*! \internal \file
   \brief
   Tests for class shared_string and helper classes

   From the stringutil tests as template.
   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   \author R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \ingroup module_data_structures
 */

#if __cplusplus >= 201103L
// for timing
#include <chrono>
// for testing thread safety in concurrent reading
#include <thread>
#endif

#include <stdexcept>
#include <iostream> // ::std::cout
#include <ios>      // ::std::left, ::std::rightusing std::cout;
#include <ostream>  // ::std::ostream (<<)
#include <iomanip>  // ::std::setw, ::std::setfill
#include <fstream>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

// define if boost::serialization is also to be compiled in and tested, add corresponding compiler flags -DHAVE_BOOST_SERIALIZATION and for linking
//#define HAVE_BOOST_SERIALIZATION
#include "gromacs/data_structures/shared_string.h"
#include "gromacs/utility/path.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

//#define DEBUG


namespace
{

class SharedStringTest : public ::testing::Test
{
    public:
        SharedStringTest()
        {
        }
        gmx::test::TestFileManager      fileManager_;
};

//! dummy stream as sink for unwanted debugging output, this is guaranteed to work/give a no-effect iostream by the C++ standard according to a post on stackoverflow.com
#if __cplusplus >= 201103L
std::ostream  cnull(nullptr);
#else
std::ostream  cnull(0);
#endif
//std::wostream wcnull(nullptr);
#ifndef DEBUG
//! By default, the debugpt points to the dummy stream cnull silencing debugging information. If macro DEBUG is defined, the debugpt points to cout (stdout) and debugging information will be printed.
std::ostream* debugpt = &cnull;
#else
//! By default, the debugpt points to the dummy stream cnull silencing debugging information. If macro DEBUG is defined, the debugpt points to cout (stdout) and debugging information will be printed.
std::ostream* debugpt = &std::cout;
#endif
//! used for stearing debug output to cout (stdout) or the dummy stream cnull
#define debug (*debugpt)

//! read strings from file to a vector of string container objects
template <class T>
void read_file(std::vector<T> &vec, const std::string &file)
{

    std::string       tmpstr;
    // open file
    std::ifstream     atxtf(file.c_str());
    // Could I open the file for reading?
    if (!atxtf.good())
    {
        std::cerr << "main: failed to open input file, " << file.c_str() << std::endl;
        throw(std::runtime_error("failed to open input file: " + file));
    }

    // line breaks are overread
    while (atxtf >> tmpstr)
    {
// the symbol table has an internal mutex, a lock-free variant would be nicer
// but seems currently not (easily) realizable, and probably not
// worthwhile unless one is really interested in text processing
        vec.push_back(tmpstr);
    }

    // close file
    atxtf.close();
}

/********************************************************************
 * Tests for simple string utilities
 */

TEST_F(SharedStringTest, ConstructDestruct)
{

    debug << std::endl;
    debug << "-------------------------------------------------------------------------------" << std::endl;
    debug << "Running test 1: creation and destruction of t_symtab, t_symbol and"              << std::endl;
    debug << "                shared_string objects ..."                                       << std::endl;
    debug << "-------------------------------------------------------------------------------" << std::endl;
    debug << std::endl;

    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test creation and destruction of helper objects:"        << std::endl;

    bool caught_exception = false;
    try
    {
        gmx::t_symbol test_symbol;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        std::cout << "Caught std::exception in constructing/destructing a \"t_symbol\" object." << std::endl;
        std::cout << "Exception message: " << e.what() << std::endl;
    }
    catch (...)
    {
        caught_exception = true;
        std::cout << "Caught non-standard exception in constructing/destructing a \"t_symbol\" object." << std::endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"t_symbol\" failed";

    caught_exception = false;
    try
    {
        const size_t test_bufsize = 80;
        const size_t test_max_entry_size = 1024;
        gmx::t_symtab     test_symtab(test_bufsize, test_max_entry_size);
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        std::cout << "Caught std::exception in constructing/destructing a \"t_symtab\" object." << std::endl;
        std::cout << "Exception message: " << e.what() << std::endl;
    }
    catch (...)
    {
        caught_exception = true;
        std::cout << "Caught non-standard exception in constructing/destructing a \"t_symtab\" object." << std::endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"t_symtab\" failed";

    caught_exception = false;
    try
    {
        gmx::shared_string test_shared_string1;
        gmx::shared_string test_shared_string2("ABC");
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        std::cout << "Caught std::exception in constructing/destructing a \"shared_string\" object." << std::endl;
        std::cout << "Exception message: " << e.what() << std::endl;
    }
    catch (...)
    {
        caught_exception = true;
        std::cout << "Caught non-standard exception in constructing/destructing a \"shared_string\" object." << std::endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"shared_string\" failed";
    debug << "-----------------------------------------------------<<" << std::endl;


    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test creation from shared_string in different variants:" << std::endl;

    gmx::shared_string test_shared_string1("");
    debug << "shared_string()              test_shared_string1 = \""   << test_shared_string1 << "\"" << std::endl;

    gmx::shared_string test_shared_string2("TestString2 A");
    debug << "shared_string(\"string\")      test_shared_string2 = \"" << test_shared_string2 << "\"" << std::endl;

    gmx::shared_string test_shared_string3(std::string("TestString3 A"));
    debug << "shared_string(std::string)   test_shared_string3 = \""   << test_shared_string3 << "\"" << std::endl;

    gmx::shared_string test_shared_string4(std::string("TestString4 A").c_str());
    debug << "shared_string(const char*)   test_shared_string4 = \""   << test_shared_string4 << "\"" << std::endl;

    gmx::shared_string test_shared_string5(gmx::shared_string("TestString5 A"));
    debug << "shared_string(shared_string) test_shared_string5 = \""   << test_shared_string5 << "\"" << std::endl;

    EXPECT_TRUE(test_shared_string1.std_str() == "" &&
                test_shared_string2.std_str() == "TestString2 A" &&
                test_shared_string3.std_str() == "TestString3 A" &&
                test_shared_string4.std_str() == "TestString4 A" &&
                test_shared_string5.std_str() == "TestString5 A")
    << "Test creation of shared string failed.";

    debug << "-----------------------------------------------------<<" << std::endl;
}

TEST_F(SharedStringTest, Assign)
{
    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test assignment of shared_string in different variants:" << std::endl;

    gmx::shared_string test_shared_string1("");
    debug << "shared_string()              test_shared_string1 = \""      << test_shared_string1 << "\"" << std::endl;

    gmx::shared_string test_shared_string2("TestString2 A");
    debug << "shared_string(\"string\")      test_shared_string2 = \""    << test_shared_string2 << "\"" << std::endl;

    gmx::shared_string test_shared_string3(std::string("TestString3 A"));
    debug << "shared_string(std::string)   test_shared_string3 = \""      << test_shared_string3 << "\"" << std::endl;

    gmx::shared_string test_shared_string4(std::string("TestString4 A").c_str());
    debug << "shared_string(const char*)   test_shared_string4 = \""      << test_shared_string4 << "\"" << std::endl;

    gmx::shared_string test_shared_string5(gmx::shared_string("TestString5 A"));
    debug << "shared_string(shared_string) test_shared_string5 = \""      << test_shared_string5 << "\"" << std::endl;

    test_shared_string1 = "TestString1 B";
    debug << "shared_string = \"string\"        test_shared_string1 = \"" << test_shared_string1 << "\"" << std::endl;

    test_shared_string2 = std::string("TestString2 B").c_str();
    debug << "shared_string = const char*     test_shared_string2 = \""   << test_shared_string2 << "\"" << std::endl;

    test_shared_string3 = std::string("TestString3 B");
    debug << "shared_string = std:string      test_shared_string3 = \""   << test_shared_string3 << "\"" << std::endl;

    test_shared_string4 = gmx::shared_string("TestString4 B");
    debug << "shared_string = shared_string   test_shared_string4 = \""   << test_shared_string4 << "\"" << std::endl;

    test_shared_string5 = gmx::shared_string("TestString5 B");
    debug << "shared_string = shared_string   test_shared_string5 = \""   << test_shared_string5 << "\"" << std::endl;

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString1 B" &&
                test_shared_string2.std_str() == "TestString2 B" &&
                test_shared_string3.std_str() == "TestString3 B" &&
                test_shared_string4.std_str() == "TestString4 B" &&
                test_shared_string5.std_str() == "TestString5 B")
    << "Test assignment of shared_string in different variants failed.";

    debug << "-----------------------------------------------------<<" << std::endl;

    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test (re)assignment to a previously existing string or another shared_string :" << std::endl;
    debug << "test_shared_string1,2,3 should share \"TestString1 D\"" << std::endl;
    debug << "test_shared_string4,5   should share \"TestString4 D\"" << std::endl;

    test_shared_string1 = "TestString1 D";
    test_shared_string2 = "TestString1 D";
    test_shared_string3 = " TestString1 D  ";
    test_shared_string4 = "TestString4 D";
    test_shared_string5 = test_shared_string4;

    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", use_count_p = " << test_shared_string1.get_use_count_p() << ", use_count_pp = " << test_shared_string1.get_use_count_pp() << ", use_count = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", use_count_p = " << test_shared_string2.get_use_count_p() << ", use_count_pp = " << test_shared_string2.get_use_count_pp() << ", use_count = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", use_count_p = " << test_shared_string3.get_use_count_p() << ", use_count_pp = " << test_shared_string3.get_use_count_pp() << ", use_count = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", use_count_p = " << test_shared_string4.get_use_count_p() << ", use_count_pp = " << test_shared_string4.get_use_count_pp() << ", use_count = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", use_count_p = " << test_shared_string5.get_use_count_p() << ", use_count_pp = " << test_shared_string5.get_use_count_pp() << ", use_count = " << test_shared_string5.get_use_count() << std::endl;

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString1 D" &&
                test_shared_string2.std_str() == "TestString1 D" &&
                test_shared_string3.std_str() == "TestString1 D" &&
                test_shared_string4.std_str() == "TestString4 D" &&
                test_shared_string5.std_str() == "TestString4 D" &&
                test_shared_string1 == test_shared_string2 &&
                test_shared_string1 == test_shared_string3 &&
                test_shared_string1 != test_shared_string4 &&
                test_shared_string4 == test_shared_string5)
    << "Test reassignment of shared_string failed.";

    debug << "-----------------------------------------------------<<" << std::endl;
}


TEST_F(SharedStringTest, ImplicitAndExplicitConversion)
{
    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test implictit and explicit conversion between char*, std::string and shared_string:" << std::endl;

    gmx::shared_string test_shared_string1 = "TestString1 B";
    debug << "shared_string = \"string\"        test_shared_string1 = \"" << test_shared_string1 << "\"" << std::endl;

    //! implict conversion from shared_string to std::string
    std::string test_string1(test_shared_string1);
    //! explicit conversion from shared_string to std::string
    std::string test_string2(test_shared_string1.std_str());
    //! conversion via explicit operator requires the explicit keyword of C++11
#if __cplusplus >= 201103L
    //! explicit conversion operator from shared_string to const char*
    std::string test_string3(static_cast<const char*>(test_shared_string1));
#else
    std::string test_string3(test_shared_string1.c_str());
#endif
    //! explicit std::string-like conversion from shared_string to const char*
    std::string test_string4(test_shared_string1.c_str());

    //! implicit conversion from shared_string to std:string
    std::string test_string5 = test_shared_string1;
    std::string test_string6 = "" + test_shared_string1;

    //! explicit conversion from shared_string to std:string
    std::string        test_string7 = test_shared_string1.std_str();
    std::string        test_string8 = "" + test_shared_string1.std_str();

    std::string        test_string9 = "TestString7 B";

    gmx::shared_string test_shared_string2 = test_string9;
    gmx::shared_string test_shared_string3 = "" + test_string9;

    EXPECT_TRUE(test_string1 == test_shared_string1 &&
                test_string2 == test_shared_string1 &&
                test_string3 == test_shared_string1 &&
                test_string4 == test_shared_string1 &&
                test_string5 == test_shared_string1 &&
                test_string6 == test_shared_string1 &&
                test_string7 == test_shared_string1 &&
                test_string8 == test_shared_string1 &&
                test_shared_string2 == test_string9 &&
                test_shared_string3 == test_string9)
    << "Test implicit/explicit conversion between shared_string and char*/std::string failed." << std::endl;

    debug << "-----------------------------------------------------<<" << std::endl;
}

TEST_F(SharedStringTest, SymbolTableIntegrity)
{
    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test comparison operators:" << std::endl;
    debug << "test_shared_string1,2,3 should share \"TestString1 C\"" << std::endl;
    debug << "test_shared_string4,5   should share \"TestString4 C\"" << std::endl;

    gmx::shared_string test_shared_string1 = "TestString1 B";
    gmx::shared_string test_shared_string2 = "TestString7 B";
    gmx::shared_string test_shared_string3 = " TestString7 B  ";
    gmx::shared_string test_shared_string4 = "TestString4 B";
    gmx::shared_string test_shared_string5 = test_shared_string4;

    test_shared_string1 = "TestString1 C";
    test_shared_string2 = "TestString1 C";
    test_shared_string3 = " TestString1 C  ";
    test_shared_string4 = "TestString4 C";
    test_shared_string5 = test_shared_string4;

    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", use_count_p = " << test_shared_string1.get_use_count_p() << ", use_count_pp = " << test_shared_string1.get_use_count_pp() << ", use_count = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", use_count_p = " << test_shared_string2.get_use_count_p() << ", use_count_pp = " << test_shared_string2.get_use_count_pp() << ", use_count = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", use_count_p = " << test_shared_string3.get_use_count_p() << ", use_count_pp = " << test_shared_string3.get_use_count_pp() << ", use_count = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", use_count_p = " << test_shared_string4.get_use_count_p() << ", use_count_pp = " << test_shared_string4.get_use_count_pp() << ", use_count = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", use_count_p = " << test_shared_string5.get_use_count_p() << ", use_count_pp = " << test_shared_string5.get_use_count_pp() << ", use_count = " << test_shared_string5.get_use_count() << std::endl;

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString1 C" &&
                test_shared_string2.std_str() == "TestString1 C" &&
                test_shared_string3.std_str() == "TestString1 C" &&
                test_shared_string4.std_str() == "TestString4 C" &&
                test_shared_string5.std_str() == "TestString4 C" &&
                test_shared_string1 == test_shared_string2 &&
                test_shared_string1 == test_shared_string3 &&
                test_shared_string1 != test_shared_string4 &&
                test_shared_string4 == test_shared_string5)
    << "Test symbol table integrity failed. Shared_string(s) have incorrect content." << std::endl;

    /*! use_count: total number of shared_string instances referencing the same symbol / string content
        use_count_p: total number of shared_pointer instances from the same set as
                     the shared_ptr of this shared_string referencing the same symbol,
                     can differ from use_count after renaming of a different string to
                     another string that also already existed in the symbol table
        use_count_pp: total number of inner shared_pointers referencing the same symbol
     */

    debug << "test_shared_string1.get_use_count() = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2.get_use_count() = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3.get_use_count() = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4.get_use_count() = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5.get_use_count() = " << test_shared_string5.get_use_count() << std::endl;
    debug << "test_shared_string1.get_use_count_p() = " << test_shared_string1.get_use_count_p() << std::endl;
    debug << "test_shared_string2.get_use_count_p() = " << test_shared_string2.get_use_count_p() << std::endl;
    debug << "test_shared_string3.get_use_count_p() = " << test_shared_string3.get_use_count_p() << std::endl;
    debug << "test_shared_string4.get_use_count_p() = " << test_shared_string4.get_use_count_p() << std::endl;
    debug << "test_shared_string5.get_use_count_p() = " << test_shared_string5.get_use_count_p() << std::endl;
    debug << "test_shared_string1.get_use_count_pp() = " << test_shared_string1.get_use_count_pp() << std::endl;
    debug << "test_shared_string2.get_use_count_pp() = " << test_shared_string2.get_use_count_pp() << std::endl;
    debug << "test_shared_string3.get_use_count_pp() = " << test_shared_string3.get_use_count_pp() << std::endl;
    debug << "test_shared_string4.get_use_count_pp() = " << test_shared_string4.get_use_count_pp() << std::endl;
    debug << "test_shared_string5.get_use_count_pp() = " << test_shared_string5.get_use_count_pp() << std::endl;

    // The C++98 version has to use a fully qualified shared_ptr also in the symbol table,
    // which is consequently also counted as instance by shared_string::get_use_count_p().
#ifdef gmx_symtab_cpp11_h
    EXPECT_TRUE(test_shared_string1.get_use_count() == test_shared_string2.get_use_count() &&
                test_shared_string1.get_use_count() == test_shared_string3.get_use_count() &&
                test_shared_string2.get_use_count() == test_shared_string3.get_use_count() &&
                test_shared_string4.get_use_count() == test_shared_string5.get_use_count() &&
                test_shared_string1.get_use_count() == test_shared_string1.get_use_count_p() &&
                test_shared_string4.get_use_count() == test_shared_string4.get_use_count_p() &&
                test_shared_string1.get_use_count() == 3  &&
                test_shared_string4.get_use_count() == 2  &&
                test_shared_string1.get_use_count_pp() == 1 &&
                test_shared_string2.get_use_count_pp() == 1 &&
                test_shared_string3.get_use_count_pp() == 1 &&
                test_shared_string4.get_use_count_pp() == 1 &&
                test_shared_string5.get_use_count_pp() == 1)
    << "Test symbol table integrity failed. Shared_string(s) have incorrect use counts. Symbol table appears corrupted." << std::endl;
#else
    EXPECT_TRUE(test_shared_string1.get_use_count() == test_shared_string2.get_use_count() &&
                test_shared_string1.get_use_count() == test_shared_string3.get_use_count() &&
                test_shared_string2.get_use_count() == test_shared_string3.get_use_count() &&
                test_shared_string4.get_use_count() == test_shared_string5.get_use_count() &&
                test_shared_string1.get_use_count() == (test_shared_string1.get_use_count_p() - 1) &&
                test_shared_string4.get_use_count() == (test_shared_string4.get_use_count_p() - 1) &&
                test_shared_string1.get_use_count() == 3  &&
                test_shared_string4.get_use_count() == 2  &&
                test_shared_string1.get_use_count_pp() == 1 &&
                test_shared_string2.get_use_count_pp() == 1 &&
                test_shared_string3.get_use_count_pp() == 1 &&
                test_shared_string4.get_use_count_pp() == 1 &&
                test_shared_string5.get_use_count_pp() == 1)
    << "Test symbol table integrity failed. Shared_string(s) have incorrect use counts. Symbol table appears corrupted." << std::endl;
#endif

    debug << "-----------------------------------------------------<<" << std::endl;
}


TEST_F(SharedStringTest, ComparisonOperators)
{
    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test comparison operators:" << std::endl;
    debug << "test_shared_string1,2,3 should share \"TestString1 C\"" << std::endl;
    debug << "test_shared_string4,5   should share \"TestString4 C\"" << std::endl;

    gmx::shared_string test_shared_string1 = "TestString1 C";
    gmx::shared_string test_shared_string2 = "TestString1 C";
    gmx::shared_string test_shared_string3 = " TestString1 C  ";
    gmx::shared_string test_shared_string4 = "TestString4 C";
    gmx::shared_string test_shared_string5 = test_shared_string4;

    // begin test of operator == ...
    EXPECT_TRUE(test_shared_string1 == test_shared_string2)
    << "Test of shared_string           == shared_string             failed."   << std::endl;
    EXPECT_TRUE(test_shared_string1.std_str() == "TestString1 C")
    << "Test of shared_string.std_str() == \"string\"                  failed." << std::endl;
    EXPECT_TRUE("TestString1 C" == test_shared_string1.std_str())
    << "Test of \"string\"                == shared_string.std_str()   failed." << std::endl;
    EXPECT_TRUE(test_shared_string1 == std::string("TestString1 C"))
    << "Test of shared_string           == std::string               failed."   << std::endl;
    EXPECT_TRUE(std::string("TestString1 C") == test_shared_string1)
    << "Test of std::string             == shared_string             failed."   << std::endl;
    EXPECT_TRUE(test_shared_string1 == "TestString1 C")
    << "Test of shared_string           == \"string\"                  failed." << std::endl;
    EXPECT_TRUE("TestString1 C" == test_shared_string1)
    << "Test of \"string\"                == shared_string             failed." << std::endl;
    // .. end test of operator ==

    // begin test of operator != ...
    EXPECT_TRUE(!(test_shared_string1 != test_shared_string2))
    << "Test of shared_string           != shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string1.std_str() != "TestString1 C"))
    << "Test of shared_string.std_str() != \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString1 C" != test_shared_string1.std_str()))
    << "Test of \"string\"                != shared_string.std_str()   failed." << std::endl;
    EXPECT_TRUE(!(test_shared_string1 != std::string("TestString1 C")))
    << "Test of shared_string           != std::string               failed."   << std::endl;
    EXPECT_TRUE(!(std::string("TestString1 C") != test_shared_string1))
    << "Test of std::string             != shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string1 != "TestString1 C"))
    << "Test of shared_string           != \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString1 C" != test_shared_string1))
    << "Test of \"string\"                != shared_string             failed." << std::endl;
    // ... end test of operator !=

    // begin test of operator > ...
    EXPECT_TRUE(!(test_shared_string1 > test_shared_string4) && (test_shared_string4 > test_shared_string1))
    << "Test of shared_string           >  shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string1.std_str() > "TestString1 C") && (test_shared_string4.std_str() > "TestString1 C"))
    << "Test of shared_string.std_str() >  \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString1 C" > test_shared_string1.std_str()) && ("TestString4 C" > test_shared_string1.std_str()))
    << "Test of \"string\"                >  shared_string.std_str()   failed." << std::endl;
    if (!(test_shared_string1 > std::string("TestString1 C")) && (test_shared_string4 > std::string("TestString1 C")))
    {
        debug << "Test of shared_string           >  std::string               passed." << std::endl;
    }
    else
    {
        debug << "Test of shared_string           >  std::string               failed." << std::endl;
    }
    EXPECT_TRUE(!(std::string("TestString1 C") > test_shared_string1) && (std::string("TestString4 C") > test_shared_string1))
    << "Test of std::string             >  shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string1 > "TestString1 C") && (test_shared_string4 > "TestString1 C"))
    << "Test of shared_string           >  \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString1 C" > test_shared_string1) && ("TestString4 C" > test_shared_string1))
    << "Test of \"string\"                >  shared_string             failed." << std::endl;
    // .. end test of operator >

    // begin test of operator < ...
    EXPECT_TRUE(!(test_shared_string4 < test_shared_string1) && (test_shared_string1 < test_shared_string4))
    << "Test of shared_string           <  shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string4.std_str() < "TestString1 C") && (test_shared_string1.std_str() < "TestString4 C"))
    << "Test of shared_string.std_str() <  \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString4 C" < test_shared_string1.std_str()) && ("TestString1 C" < test_shared_string4.std_str()))
    << "Test of \"string\"                <  shared_string.std_str()   failed." << std::endl;
    if (!(test_shared_string4 < std::string("TestString1 C")) && (test_shared_string1 < std::string("TestString4 C")))
    {
        debug << "Test of shared_string           <  std::string               passed." << std::endl;
    }
    else
    {
        debug << "Test of shared_string           <  std::string               failed." << std::endl;
    }
    EXPECT_TRUE(!(std::string("TestString4 C") < test_shared_string1) && (std::string("TestString1 C") < test_shared_string4))
    << "Test of std::string             <  shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string4 < "TestString1 C") && (test_shared_string1 < "TestString4 C"))
    << "Test of shared_string           <  \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString4 C" < test_shared_string1) && ("TestString1 C" < test_shared_string4))
    << "Test of \"string\"                <  shared_string             failed." << std::endl;
    // .. end test of operator <

    // begin test of operator >= ...
    EXPECT_TRUE(!(test_shared_string1 >= test_shared_string4) && (test_shared_string4 >= test_shared_string1) && (test_shared_string2 >= test_shared_string1))
    << "Test of shared_string           >= shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string1.std_str() >= "TestString4 C") && (test_shared_string4.std_str() >= "TestString1 C") && (test_shared_string2.std_str() >= "TestString1 C"))
    << "Test of shared_string.std_str() >= \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString1 C" >= test_shared_string4.std_str()) && ("TestString4 C" >= test_shared_string1.std_str()) && ("TestString1 C" >= test_shared_string1.std_str()))
    << "Test of \"string\"                >= shared_string.std_str()   failed." << std::endl;
    EXPECT_TRUE(!(test_shared_string1 >= std::string("TestString4 C")) && (test_shared_string4 >= std::string("TestString1 C")) && (test_shared_string2 >= std::string("TestString1 C")))
    << "Test of shared_string           >= std::string               failed."   << std::endl;
    EXPECT_TRUE(!(std::string("TestString1 C") >= test_shared_string4) && (std::string("TestString4 C") >= test_shared_string1) && (std::string("TestString1 C") >= test_shared_string1))
    << "Test of std::string             >= shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string1 >= "TestString4 C") && (test_shared_string4 >= "TestString1 C") && (test_shared_string2 >= "TestString1 C"))
    << "Test of shared_string           >= \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString1 C" >= test_shared_string4) && ("TestString4 C" >= test_shared_string1) && ("TestString1 C" >= test_shared_string1))
    << "Test of \"string\"                >= shared_string             failed." << std::endl;
    // .. end test of operator >=

    // begin test of operator <= ...
    EXPECT_TRUE(!(test_shared_string4 <= test_shared_string1) && (test_shared_string1 <= test_shared_string4) && (test_shared_string1 <= test_shared_string2))
    << "Test of shared_string           <= shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string4.std_str() <= "TestString1 C") && (test_shared_string1.std_str() <= "TestString4 C") && (test_shared_string1.std_str() <= "TestString1 C"))
    << "Test of shared_string.std_str() <= \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString4 C" <= test_shared_string1.std_str()) && ("TestString1 C" <= test_shared_string4.std_str()) && ("TestString1 C" <= test_shared_string2.std_str()))
    << "Test of \"string\"                <= shared_string.std_str()   failed." << std::endl;
    EXPECT_TRUE(!(test_shared_string4 <= std::string("TestString1 C")) && (test_shared_string1 <= std::string("TestString4 C")) && (test_shared_string1 <= std::string("TestString1 C")))
    << "Test of shared_string           <= std::string               failed."   << std::endl;
    EXPECT_TRUE(!(std::string("TestString4 C") <= test_shared_string1) && (std::string("TestString1 C") <= test_shared_string4) && (std::string("TestString1 C") <= test_shared_string2))
    << "Test of std::string             <= shared_string             failed."   << std::endl;
    EXPECT_TRUE(!(test_shared_string4 <= "TestString1 C") && (test_shared_string1 <= "TestString4 C") && (test_shared_string1 <= "TestString1 C"))
    << "Test of shared_string           <= \"string\"                  failed." << std::endl;
    EXPECT_TRUE(!("TestString4 C" <= test_shared_string1) && ("TestString1 C" <= test_shared_string4) && ("TestString1 C" <= test_shared_string1))
    << "Test of \"string\"                <= shared_string             failed." << std::endl;
    // .. end test of operator <=

    debug << "-----------------------------------------------------<<" << std::endl;
}

TEST_F(SharedStringTest, ReplaceAll)
{
    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test renaming of a previously assigned string via shared_string member function with single argument:" << std::endl;
    debug << "test_shared_string1,2,3,4,5 should share \"TestString4 E\"" << std::endl;

    gmx::shared_string test_shared_string1("TestString1 E");
    gmx::shared_string test_shared_string2("TestString1 E");
    gmx::shared_string test_shared_string3(" TestString1 E  ");
    gmx::shared_string test_shared_string4("TestString4 E");
    gmx::shared_string test_shared_string5(test_shared_string4);

    debug << "initial strings before renaming:" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", use_count_p = " << test_shared_string1.get_use_count_p() << ", use_count_pp = " << test_shared_string1.get_use_count_pp() << ", use_count = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", use_count_p = " << test_shared_string2.get_use_count_p() << ", use_count_pp = " << test_shared_string2.get_use_count_pp() << ", use_count = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", use_count_p = " << test_shared_string3.get_use_count_p() << ", use_count_pp = " << test_shared_string3.get_use_count_pp() << ", use_count = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", use_count_p = " << test_shared_string4.get_use_count_p() << ", use_count_pp = " << test_shared_string4.get_use_count_pp() << ", use_count = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", use_count_p = " << test_shared_string5.get_use_count_p() << ", use_count_pp = " << test_shared_string5.get_use_count_pp() << ", use_count = " << test_shared_string5.get_use_count() << std::endl;
    if (debugpt != &cnull)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debug_stats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debug_stats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debug_stats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debug_stats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debug_stats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debug_symtab_stats();
    }

    test_shared_string4.replace_all("TestString1 E");

    debug << "final strings after renaming via test_shared_string4.replace_all(\"TestString1 E\"):" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", use_count_p = " << test_shared_string1.get_use_count_p() << ", use_count_pp = " << test_shared_string1.get_use_count_pp() << ", use_count = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", use_count_p = " << test_shared_string2.get_use_count_p() << ", use_count_pp = " << test_shared_string2.get_use_count_pp() << ", use_count = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", use_count_p = " << test_shared_string3.get_use_count_p() << ", use_count_pp = " << test_shared_string3.get_use_count_pp() << ", use_count = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", use_count_p = " << test_shared_string4.get_use_count_p() << ", use_count_pp = " << test_shared_string4.get_use_count_pp() << ", use_count = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", use_count_p = " << test_shared_string5.get_use_count_p() << ", use_count_pp = " << test_shared_string5.get_use_count_pp() << ", use_count = " << test_shared_string5.get_use_count() << std::endl;
    if (debugpt != &cnull)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debug_stats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debug_stats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debug_stats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debug_stats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debug_stats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debug_symtab_stats();
    }

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString1 E" &&
                test_shared_string2.std_str() == "TestString1 E" &&
                test_shared_string3.std_str() == "TestString1 E" &&
                test_shared_string4.std_str() == "TestString1 E" &&
                test_shared_string5.std_str() == "TestString1 E" &&
                test_shared_string1 == test_shared_string2 &&
                test_shared_string1 == test_shared_string3 &&
                test_shared_string1 == test_shared_string4 &&
                test_shared_string1 == test_shared_string5)
    << "Test replace_all() with target \"string\" as only argument failed. Shared_string(s) have incorrect content." << std::endl;

    /*! use_count: total number of shared_string instances referencing the same symbol / string content
        use_count_p: total number of shared_pointer instances from the same set as
                     the shared_ptr of this shared_string referencing the same symbol,
                     can differ from use_count after renaming of a different string to
                     another string that also already existed in the symbol table
        use_count_pp: total number of inner shared_pointers referencing the same symbol
     */

    EXPECT_TRUE(test_shared_string1.get_use_count() == 5)
    << "Test replace_all() with target \"string\" as only argument failed. Shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;

    debug << "-----------------------------------------------------<<" << std::endl;


    debug << ">>-----------------------------------------------------" << std::endl;
    debug << "Test renaming of a previously assigned string via free function with \"before\" and \"after\" arguments:" << std::endl;
    debug << "test_shared_string1,2,3,4,5 should share \"TestString4 E\"" << std::endl;

    test_shared_string1 = "TestString1 E";
    test_shared_string2 = "TestString1 E";
    test_shared_string3 = " TestString1 E  ";
    test_shared_string4 = "TestString4 E";
    test_shared_string5 = test_shared_string4;

    debug << "initial strings before renaming:" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", use_count_p = " << test_shared_string1.get_use_count_p() << ", use_count_pp = " << test_shared_string1.get_use_count_pp() << ", use_count = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", use_count_p = " << test_shared_string2.get_use_count_p() << ", use_count_pp = " << test_shared_string2.get_use_count_pp() << ", use_count = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", use_count_p = " << test_shared_string3.get_use_count_p() << ", use_count_pp = " << test_shared_string3.get_use_count_pp() << ", use_count = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", use_count_p = " << test_shared_string4.get_use_count_p() << ", use_count_pp = " << test_shared_string4.get_use_count_pp() << ", use_count = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", use_count_p = " << test_shared_string5.get_use_count_p() << ", use_count_pp = " << test_shared_string5.get_use_count_pp() << ", use_count = " << test_shared_string5.get_use_count() << std::endl;
    if (debugpt != &cnull)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debug_stats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debug_stats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debug_stats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debug_stats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debug_stats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debug_symtab_stats();
    }

    gmx::replace_all_shared_strings("TestString1 E", "TestString4 E");

    debug << "final strings after renaming via replace_all_shared_strings(\"TestString1 E\", \"TestString4 E\")" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", use_count_p = " << test_shared_string1.get_use_count_p() << ", use_count_pp = " << test_shared_string1.get_use_count_pp() << ", use_count = " << test_shared_string1.get_use_count() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", use_count_p = " << test_shared_string2.get_use_count_p() << ", use_count_pp = " << test_shared_string2.get_use_count_pp() << ", use_count = " << test_shared_string2.get_use_count() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", use_count_p = " << test_shared_string3.get_use_count_p() << ", use_count_pp = " << test_shared_string3.get_use_count_pp() << ", use_count = " << test_shared_string3.get_use_count() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", use_count_p = " << test_shared_string4.get_use_count_p() << ", use_count_pp = " << test_shared_string4.get_use_count_pp() << ", use_count = " << test_shared_string4.get_use_count() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", use_count_p = " << test_shared_string5.get_use_count_p() << ", use_count_pp = " << test_shared_string5.get_use_count_pp() << ", use_count = " << test_shared_string5.get_use_count() << std::endl;
    if (debugpt != &cnull)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debug_stats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debug_stats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debug_stats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debug_stats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debug_stats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debug_symtab_stats();
    }

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString4 E" &&
                test_shared_string2.std_str() == "TestString4 E" &&
                test_shared_string3.std_str() == "TestString4 E" &&
                test_shared_string4.std_str() == "TestString4 E" &&
                test_shared_string5.std_str() == "TestString4 E" &&
                test_shared_string1 == test_shared_string2 &&
                test_shared_string1 == test_shared_string3 &&
                test_shared_string1 == test_shared_string4 &&
                test_shared_string1 == test_shared_string5)
    << "Test replace_all() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect content." << std::endl;

    /*! use_count: total number of shared_string instances referencing the same symbol / string content
        use_count_p: total number of shared_pointer instances from the same set as
                     the shared_ptr of this shared_string referencing the same symbol,
                     can differ from use_count after renaming of a different string to
                     another string that also already existed in the symbol table
        use_count_pp: total number of inner shared_pointers referencing the same symbol
     */

    EXPECT_TRUE(test_shared_string1.get_use_count() == 5)
    << "Test replace_all() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;

    debug << "-----------------------------------------------------<<" << std::endl;
}

TEST_F(SharedStringTest, ReadFromFileSingleThreaded)
{
    debug << ">>--------------------------------------------------------------" << std::endl;
    debug << "Test single-threaded reading from text file into shared_strings " << std::endl;
    debug << "and test some functions (insert/retrieve/remove/rename) entries," << std::endl;
    debug << "report (approximate) total amount of occupied memory) ..."        << std::endl;

    std::string textfile("refdata/textfile.1.txt");
    if (textfile != "")
    {

        std::vector<std::string>        string_vec;
        std::vector<gmx::shared_string> shared_vec;
        std::string                     tmpstr;
        size_t                          tmp_nstrings = 0;

        // open file
        std::ifstream atxtf(fileManager_.getInputFilePath(textfile.c_str()).c_str());
        // Could I open the file for reading?
        EXPECT_TRUE(atxtf.good())
        << "Failed to open file \"" << textfile << "\" for reading" << std::endl;
        if (!atxtf.good())
        {
            std::cerr << "SharedStringTest ReadFromFileSingleThreaded: failed to open input file, " << textfile.c_str() << std::endl;
            return;
        }

        // linebreaks are overread
        while (atxtf >> tmpstr)
        {
            tmp_nstrings++;
            string_vec.push_back(tmpstr);
            shared_vec.push_back(tmpstr);
        }

        // close file
        atxtf.close();

        EXPECT_TRUE(string_vec.size() == shared_vec.size()) << "Single-threaded reading failed, differing numbers of strings and shared_strings" << std::endl;

        bool flag = true;
        for (size_t i = 0; i < string_vec.size(); ++i)
        {
            if (string_vec[i] != shared_vec[i])
            {
                debug << "Comparison of contents read to string and shared_string failed:" << std::endl;
                debug << "shared_vec[" << i << "] = " << shared_vec[i] << " while " << std::endl
                << "string_vec[" << i << "] = " << string_vec[i] << std::endl;
                flag = false;
                break;
            }
        }
        EXPECT_TRUE(flag) << "Single-threaded reading failed, differing numbers or contents of strings and shared_strings" << std::endl;

        debug << "Number of string(s) read = " << tmp_nstrings << std::endl;
        if (shared_vec.size() != 0)
        {
            debug << "Number of unique string(s) read = " << shared_vec[shared_vec.size()-1].get_nr() << std::endl;
            debug << "Size of shared_string infrastructure = " << sizeof(gmx::shared_string) << std::endl;
            debug << "Memory occupied by symbol table (estimated) = " << shared_vec[0].get_total_mem() << " byte" << std::endl;
        }

        if (shared_vec.size() != 0)
        {
            std::string tmpstr1 = "What_A_Terribly_Looooong_Name_For_A_Binding_Form";
            std::string tmpstr2("A_Shorter_Name");
//            std::string tmpstr2(tmpstr1);
            // remove preexisting instances of the strings to stick to the simple case of one instance only
            gmx::replace_all_shared_strings(tmpstr1, "ABC");
            gmx::replace_all_shared_strings(tmpstr2, "ABC");
            shared_vec.push_back(gmx::shared_string("What_A_Terribly_Looooong_Name_For_A_Binding_Form"));

            shared_vec[shared_vec.size()-1] = tmpstr1;
            ssize_t mem_before = shared_vec[shared_vec.size()-1].get_total_mem();
            gmx::replace_all_shared_strings(tmpstr1, tmpstr2);
            // renaming should have been tested already
            if (shared_vec[shared_vec.size()-1] == tmpstr2)
            {
                debug << "Successfully renamed entry \"" << tmpstr1 << "\" to \"" << tmpstr2 << "\" in the symbol table" << std::endl;
            }
            else
            {
                debug << "Could not rename entry \"" << tmpstr1 << "\" to \"" << tmpstr2  << "\" in the symbol table." << std::endl;
                debug << " Instead, got: " << shared_vec[shared_vec.size()-1] << std::endl;
            }

            ssize_t mem_after = shared_vec[shared_vec.size()-1].get_total_mem();
            EXPECT_TRUE( ((static_cast<ssize_t>(tmpstr1.length()) - static_cast<ssize_t>(tmpstr2.length())) * static_cast<ssize_t>(sizeof(char)))
                       >= (static_cast<ssize_t>(mem_after)        - static_cast<ssize_t>(mem_before)) )
            << "Memory book keeping upon renaming seems not to work. difference expected = "
            << ((static_cast<ssize_t>(tmpstr2.length()) - static_cast<ssize_t>(tmpstr1.length())) * static_cast<ssize_t>(sizeof(char)))
            << " got = " << (static_cast<ssize_t>(mem_after) - static_cast<ssize_t>(mem_before))
            << std::endl;

            // test symbol table cleanup
            // introduce two unique strings at the end, and remove the first one
            // to create a vacant position / gap in the symbol table for clean()
            shared_vec.push_back("UniqueTestSharedStringXYZxyzABCabc-0");
            shared_vec.push_back("UniqueTestSharedStringXYZxyzABCabc-1");
            shared_vec[shared_vec.size()-2] = shared_vec[0];
            mem_before = shared_vec[0].get_total_mem();
            shared_vec[0].clean();
            mem_after = shared_vec[0].get_total_mem();
            EXPECT_TRUE((static_cast<ssize_t>(mem_after) - static_cast<ssize_t>(mem_before)) <= 0)
            << "Memory book keeping and/or cleanup() seems not to work."
            << "Expected reduction of memory foortprint by >= 0 bytes,"
            << " got = " << (static_cast<ssize_t>(mem_after) - static_cast<ssize_t>(mem_before))
            << std::endl;


            // after the cleanup, we should still be in possession of all strings
            bool flag = true;
            for (size_t i = 0; i < string_vec.size(); ++i)
            {
                if (string_vec[i] != shared_vec[i])
                {
                    debug << "Comparison of string and shared_string after cleaning up the symbol table failed:" << std::endl;
                    debug << "shared_vec[" << i << "] = " << shared_vec[i] << " while " << std::endl
                    << "string_vec[" << i << "] = " << string_vec[i] << std::endl;
                    flag = false;
                    break;
                }
            }
            EXPECT_TRUE(flag) << "Cleanup after single-threaded reading failed, differing numbers or contents of strings and shared_strings." << std::endl;

        }
    }

    debug << "-----------------------------------------------------<<" << std::endl;
}

// the C++98 version is also usable with multi-threading and passes this test,
// the distinction is made here to check for availability of C++11 std::threads.
#if __cplusplus >= 201103L

//! this one is intended for testing with helgrind or similar for threading errors
//! with GNU libstdc++, setting the environment variable GLIBCPP_FORCE_NEW to 1
//! avoids false-positives in detecting data races due to reallocating memory in
//! regions previously allocated and freed again by other threads
TEST_F(SharedStringTest, ReadFromFileMultiThreaded)
{
    debug << ">>--------------------------------------------------------------" << std::endl;
    debug << "Test multi-threaded reading from text file into shared_strings  " << std::endl;
    debug << "and test some functions (insert/retrieve/remove/rename) entries," << std::endl;
    debug << "report (approximate) total amount of occupied memory) ..."        << std::endl;

    const size_t nthreads = (std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 2);
    std::vector<std::thread> tvec;
    debug << "Number of threads used = " << nthreads << std::endl;

    const size_t ntextfiles = 3;
    debug << "Number of textfiles assigned round-robin to the reading threads = " << ntextfiles << std::endl;
    std::vector<std::string> textfiles;
    for (size_t i = 0; i < ntextfiles; ++i)
    {
        char buffer[50];
        sprintf(buffer, "refdata/textfile.%li.txt", i + 1);
        textfiles.push_back(std::string(buffer));
    }
    // map ntextfiles input files to nthreads threads for tests with arbitrary numbers of threads,
    // regardless of the number of available text files
    bool all_resolved = true;
    std::vector<std::string> tfs_resolved;
    for (size_t i = 0; i < nthreads; ++i)
    {
        std::string tmpstr1(textfiles[(i % ntextfiles)]);
        std::string tmpstr2(fileManager_.getInputFilePath(tmpstr1.c_str()));
        if (tmpstr2 != "")
        {
            tfs_resolved.push_back(tmpstr2);
        }
        else
        {
            all_resolved = false;
            break;
        }
    }

    if (all_resolved)
    {

        std::vector<std::string>                      string_vec;
        std::vector<std::vector<gmx::shared_string> > shared_vec_threads;
        std::vector<std::vector<std::string> >        string_vec_threads;
        shared_vec_threads.resize(nthreads);
        string_vec_threads.resize(nthreads);
        std::vector<gmx::shared_string>               shared_vec;
        std::string tmpstr;

        for (size_t i = 0; i < nthreads; ++i)
        {
            tvec.push_back( std::thread(read_file<gmx::shared_string>, std::ref(shared_vec_threads[i]), std::ref(tfs_resolved[i])) );
        }
        for (size_t i = 0; i < nthreads; ++i)
        {
            tvec[i].join();
        }
        tvec.clear();

        size_t n_shared_strings = 0;
        for (std::vector<std::vector<gmx::shared_string> >::const_iterator it = shared_vec_threads.begin(); it != shared_vec_threads.end(); ++it)
        {
            n_shared_strings += it->size();
        }
        // collect all shared strings
        shared_vec.reserve( n_shared_strings ); // preallocate memory
        for (std::vector<std::vector<gmx::shared_string> >::iterator it = shared_vec_threads.begin(); it != shared_vec_threads.end(); ++it)
        {
            shared_vec.insert(shared_vec.end(), it->begin(), it->end());
            it->clear();
        }

        for (size_t i = 0; i < nthreads; ++i)
        {
            tvec.push_back( std::thread(read_file<std::string>, std::ref(string_vec_threads[i]), std::ref(tfs_resolved[i])) );
        }
        for (size_t i = 0; i < nthreads; ++i)
        {
            tvec[i].join();
        }
        tvec.clear();

        size_t n_std_strings = 0;
        for (std::vector<std::vector<std::string> >::const_iterator it = string_vec_threads.begin(); it != string_vec_threads.end(); ++it)
        {
            n_std_strings += it->size();
        }
        // collect all shared strings
        shared_vec.reserve( n_shared_strings ); // preallocate memory
        for (std::vector<std::vector<std::string> >::iterator it = string_vec_threads.begin(); it != string_vec_threads.end(); ++it)
        {
            string_vec.insert(string_vec.end(), it->begin(), it->end());
            it->clear();
        }

        EXPECT_TRUE(string_vec.size() == shared_vec.size()) << "Multi-threaded reading failed, differing numbers of strings and shared_strings" << std::endl;

        bool flag = true;
        for (size_t i = 0; i < string_vec.size(); ++i)
        {
            if (string_vec[i] != shared_vec[i])
            {
                debug << "Comparison of contents read to string and shared_string failed:" << std::endl;
                debug << "shared_vec[" << i << "] = " << shared_vec[i] << " while " << std::endl
                << "string_vec[" << i << "] = " << string_vec[i] << std::endl;
                flag = false;
                break;
            }
        }
        EXPECT_TRUE(flag) << "Multi-threaded reading failed, differing numbers or contents of strings and shared_strings" << std::endl;

        debug << "Number of string(s) read = " << shared_vec.size() << std::endl;
        if (shared_vec.size() != 0)
        {
            debug << "Number of unique string(s) read = " << shared_vec[shared_vec.size()-1].get_nr() << std::endl;
            debug << "Size of shared_string infrastructure = " << sizeof(gmx::shared_string) << std::endl;
            debug << "Memory occupied by symbol table (estimated) = " << shared_vec[0].get_total_mem() << " byte" << std::endl;
        }
    }

    debug << "-----------------------------------------------------<<" << std::endl;
}

#endif // ... __cplusplus >= 201103L, C++11 dependent  multithreading test

} // namespace
