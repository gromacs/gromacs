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
/*! \internal \file
   \brief
   Tests for class shared_string and helper classes

   From the stringutil tests as template.
   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   \author R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/data_structures/shared_string.h"

#include <thread>

#include "data_struct_test_commons.h"

namespace
{

// older(?) MS Visual C++ compilers don' know ssize_t
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

class SharedStringTest : public ::testing::Test
{
    public:
        SharedStringTest()
            : debug(nullptr), error(nullptr)
        {
        }
        void SetUp()
        {
            // trigger debug output if requested via commandline flag -debug
            if (gmx::data_struct_test::g_debugBool)
            {
                activateDebugOut();
                activateErrorOut();
            }
        }
        void TearDown()
        {
            deactivateDebugOut();
            deactivateErrorOut();
        }

        void activateDebugOut()
        {
            debug.rdbuf(std::cout.rdbuf());
        }
        void activateErrorOut()
        {
            error.rdbuf(std::cerr.rdbuf());
        }
        void deactivateDebugOut()
        {
            debug.rdbuf(nullptr);
        }
        void deactivateErrorOut()
        {
            error.rdbuf(nullptr);
        }

        gmx::test::TestFileManager   fileManager_;
    public:
        // debugging output stream, silent by default
        std::ostream                 debug;
        // error output stream, silent by default
        std::ostream                 error;
};


//! read strings from file to a vector of string container objects
//! \param[in]   vec    vector for inserting the read strings
//! \param[in]   file   file name
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

/*********************************************************************/


TEST_F(SharedStringTest, ConstructDestructSharedStringHelperObjects)
{
    EXPECT_NO_THROW_GMX({ gmx::Symbol test_symbol; } )
    << "Construction destruction of \"Symbol\" failed";

    EXPECT_NO_THROW_GMX({ gmx::Symbuf test_symbuf; } )
    << "Construction destruction of \"Symbuf\" failed";


    EXPECT_NO_THROW_GMX(
            {
                const size_t    test_bufsize        = 80;
                const size_t    test_max_entry_size = 1024;
                gmx::Symtab     test_symtab(test_bufsize, test_max_entry_size);
            }
            ) << "Construction destruction of \"Symtab\" failed";

    EXPECT_NO_THROW_GMX(
            {
                gmx::shared_string test_shared_string1;
                gmx::shared_string test_shared_string2("ABC");
                gmx::shared_string test_shared_string3("ABC", -1);
                gmx::shared_string test_shared_string4("ABC",  2);
            }
            ) << "Construction destruction of \"shared_string\" failed";
}

TEST_F(SharedStringTest, ConstructDestructSharedString)
{
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
}

TEST_F(SharedStringTest, Assign)
{
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

    debug << "Test (re)assignment to a previously existing string or another shared_string :" << std::endl;
    debug << "test_shared_string1,2,3 should share \"TestString1 D\"" << std::endl;
    debug << "test_shared_string4,5   should share \"TestString4 D\"" << std::endl;

    test_shared_string1 = "TestString1 D";
    test_shared_string2 = "TestString1 D";
    test_shared_string3 = " TestString1 D  ";
    test_shared_string4 = "TestString4 D";
    test_shared_string5 = test_shared_string4;

    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", useCount = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", useCount = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", useCount = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", useCount = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", useCount = " << test_shared_string5.getUseCount() << std::endl;

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
}


TEST_F(SharedStringTest, ImplicitAndExplicitConversion)
{
    debug << "Test implicit and explicit conversion between char*, std::string and shared_string:" << std::endl;

    gmx::shared_string test_shared_string1 = "TestString1 B";
    debug << "shared_string = \"string\"        test_shared_string1 = \"" << test_shared_string1 << "\"" << std::endl;

    //! implict conversion from shared_string to std::string
    std::string test_string1(test_shared_string1);
    //! explicit conversion from shared_string to std::string
    std::string test_string2(test_shared_string1.std_str());
    //! explicit conversion operator from shared_string to const char*
    std::string test_string3(static_cast<const char*>(test_shared_string1));
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

}

TEST_F(SharedStringTest, SymbolTableIntegrity)
{
    debug << "Testing symbol table integrity:"                        << std::endl;
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

    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", useCount = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", useCount = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", useCount = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", useCount = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", useCount = " << test_shared_string5.getUseCount() << std::endl;

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
        useCountP: total number of shared_pointer instances from the same set as
                     the shared_ptr of this shared_string referencing the same symbol,
                     can differ from use_count after renaming of a different string to
                     another string that also already existed in the symbol table
        useCountPP: total number of inner shared_pointers referencing the same symbol
     */

    debug << "test_shared_string1.getUseCount() = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2.getUseCount() = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3.getUseCount() = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4.getUseCount() = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5.getUseCount() = " << test_shared_string5.getUseCount() << std::endl;
    debug << "test_shared_string1.getUseCountP() = " << test_shared_string1.getUseCountP() << std::endl;
    debug << "test_shared_string2.getUseCountP() = " << test_shared_string2.getUseCountP() << std::endl;
    debug << "test_shared_string3.getUseCountP() = " << test_shared_string3.getUseCountP() << std::endl;
    debug << "test_shared_string4.getUseCountP() = " << test_shared_string4.getUseCountP() << std::endl;
    debug << "test_shared_string5.getUseCountP() = " << test_shared_string5.getUseCountP() << std::endl;
    debug << "test_shared_string1.getUseCountPP() = " << test_shared_string1.getUseCountPP() << std::endl;
    debug << "test_shared_string2.getUseCountPP() = " << test_shared_string2.getUseCountPP() << std::endl;
    debug << "test_shared_string3.getUseCountPP() = " << test_shared_string3.getUseCountPP() << std::endl;
    debug << "test_shared_string4.getUseCountPP() = " << test_shared_string4.getUseCountPP() << std::endl;
    debug << "test_shared_string5.getUseCountPP() = " << test_shared_string5.getUseCountPP() << std::endl;

    EXPECT_TRUE(test_shared_string1.getUseCount() == test_shared_string2.getUseCount() &&
                test_shared_string1.getUseCount() == test_shared_string3.getUseCount() &&
                test_shared_string2.getUseCount() == test_shared_string3.getUseCount() &&
                test_shared_string4.getUseCount() == test_shared_string5.getUseCount() &&
                test_shared_string1.getUseCount() == test_shared_string1.getUseCountP() &&
                test_shared_string4.getUseCount() == test_shared_string4.getUseCountP() &&
                test_shared_string1.getUseCount() == 3  &&
                test_shared_string4.getUseCount() == 2  &&
                test_shared_string1.getUseCountPP() == 1 &&
                test_shared_string2.getUseCountPP() == 1 &&
                test_shared_string3.getUseCountPP() == 1 &&
                test_shared_string4.getUseCountPP() == 1 &&
                test_shared_string5.getUseCountPP() == 1)
    << "Test symbol table integrity failed. Shared_string(s) have incorrect use counts. Symbol table appears corrupted." << std::endl;

}


TEST_F(SharedStringTest, ComparisonOperators)
{
    debug << "Testing comparison operators:" << std::endl;
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
}

TEST_F(SharedStringTest, ReplaceAll)
{
    debug << "Test renaming of a previously assigned string via shared_string member function with single argument:" << std::endl;
    debug << "test_shared_string1,2,3,4,5 should share \"TestString4 E\"" << std::endl;

    gmx::shared_string test_shared_string1("TestString1 E");
    gmx::shared_string test_shared_string2("TestString1 E");
    gmx::shared_string test_shared_string3(" TestString1 E  ");
    gmx::shared_string test_shared_string4("TestString4 E");
    gmx::shared_string test_shared_string5(test_shared_string4);

    debug << "initial strings before renaming:" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", useCount = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", useCount = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", useCount = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", useCount = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", useCount = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
    }

    test_shared_string4.replaceAll("TestString1 E");

    debug << "final strings after renaming via test_shared_string4.replaceAll(\"TestString1 E\"):" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", useCount = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", useCount = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", useCount = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", useCount = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", useCount = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
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
    << "Test replaceAll() with target \"string\" as only argument failed. Shared_string(s) have incorrect content." << std::endl;

    /*! useCount  : total number of shared_string instances referencing the same Symbol / string content
        useCountP : total number of shared_ptr instances from the same set as
                    the shared_ptr of this shared_string referencing the same Symbol,
                    can differ from use_count after renaming of a different string to
                    another string that also already existed in the symbol table
        useCountPP: total number of inner shared_ptrs referencing the same Symbol
     */

    EXPECT_TRUE(test_shared_string1.getUseCount() == 5)
    << "Test replaceAll() with target \"string\" as only argument failed. shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;

    debug << "Test renaming of a previously assigned string via free function with \"before\" and \"after\" arguments:" << std::endl;
    debug << "test_shared_string1,2,3,4,5 should share \"TestString4 E\"" << std::endl;

    test_shared_string1 = "TestString1 E";
    test_shared_string2 = "TestString1 E";
    test_shared_string3 = " TestString1 E  ";
    test_shared_string4 = "TestString4 E";
    test_shared_string5 = test_shared_string4;

    debug << "initial strings before renaming:" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", use_count = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", use_count = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", use_count = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", use_count = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", use_count = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
    }

    gmx::replaceAllSharedStrings("TestString1 E", "TestString4 E");

    debug << "final strings after renaming via replaceAllSharedStrings(\"TestString1 E\", \"TestString4 E\")" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", use_count = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", use_count = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", use_count = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", use_count = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", use_count = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
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
    << "Test replaceAll() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect content." << std::endl;

    /*! useCount  : total number of shared_string instances referencing the same symbol / string content
        useCountP : total number of shared_ptr instances from the same set as
                    the shared_ptr of this shared_string referencing the same symbol,
                    can differ from useCount after renaming of a different string to
                    another string that also already existed in the symbol table
        useCountPP: total number of inner shared_ptrs referencing the same symbol
     */

    EXPECT_TRUE(test_shared_string1.getUseCount() == 5)
    << "Test replaceAll() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;
}

TEST_F(SharedStringTest, ReplaceAllDistinctStringGroups)
{
    debug << "Test renaming of a previously assigned string via shared_string member function with single argument:" << std::endl;
    debug << "test_shared_string1,2,3 should share \"TestString4 E\"" << std::endl;

    gmx::shared_string test_shared_string1("TestString4 E", 1);
    gmx::shared_string test_shared_string2("TestString4 E", 1);
    gmx::shared_string test_shared_string3(" TestString4 E  ", 1);
    gmx::shared_string test_shared_string4("TestString4 E", 2);
    // should be assigned the default group ID -1
    gmx::shared_string test_shared_string5("TestString4 E");

    debug << "initial strings before renaming:" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", useCount = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", useCount = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", useCount = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", useCount = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", useCount = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
    }

    test_shared_string3.replaceAll("TestString1 E");

    debug << "final strings after renaming via test_shared_string3.replaceAll(\"TestString1 E\"):" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", useCount = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", useCount = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", useCount = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", useCount = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", useCount = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
    }

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString1 E" &&
                test_shared_string2.std_str() == "TestString1 E" &&
                test_shared_string3.std_str() == "TestString1 E" &&
                test_shared_string4.std_str() == "TestString4 E" &&
                test_shared_string5.std_str() == "TestString4 E" &&
                test_shared_string1 == test_shared_string2 &&
                test_shared_string1 == test_shared_string3 &&
                test_shared_string4 == test_shared_string5)
    << "Test replaceAll() with target \"string\" as only argument failed. Shared_string(s) have incorrect content." << std::endl;

    /*! useCount  : total number of shared_string instances referencing the same Symbol / string content
        useCountP : total number of shared_ptr instances from the same set as
                    the shared_ptr of this shared_string referencing the same Symbol,
                    can differ from use_count after renaming of a different string to
                    another string that also already existed in the symbol table
        useCountPP: total number of inner shared_ptrs referencing the same Symbol
     */

    EXPECT_TRUE(test_shared_string1.getGroup() ==  1 &&
                test_shared_string4.getGroup() ==  2 &&
                test_shared_string5.getGroup() == -1)
    << "Test replaceAll() with target \"string\" as only argument failed. shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;

    EXPECT_TRUE(test_shared_string1.getUseCount() == 3 &&
                test_shared_string4.getUseCount() == 1 &&
                test_shared_string5.getUseCount() == 1)
    << "Test replaceAll() with target \"string\" as only argument failed. shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;


    debug << "Test renaming of a previously assigned string via free function with \"before\" and \"after\" arguments:" << std::endl;
    debug << "test_shared_string1,2,3,4 should share \"TestString4 E\"" << std::endl;

    test_shared_string1 = "TestString1 E";
    test_shared_string2 = "TestString1 E";
    test_shared_string3 = " TestString1 E  ";
    test_shared_string4 = "  TestString1 E";
    test_shared_string2.setGroup(test_shared_string1.getGroup());
    test_shared_string3.setGroup(test_shared_string1.getGroup());
    test_shared_string4.setGroup(test_shared_string1.getGroup());
    test_shared_string5 = test_shared_string1;
    test_shared_string5.setGroup(test_shared_string1.getGroup() + 1);

    debug << "initial strings before renaming:" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", group = " << test_shared_string1.getGroup() << ", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", use_count = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", group = " << test_shared_string2.getGroup() << ", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", use_count = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", group = " << test_shared_string3.getGroup() << ", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", use_count = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", group = " << test_shared_string4.getGroup() << ", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", use_count = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", group = " << test_shared_string5.getGroup() << ", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", use_count = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
    }

    long int old_gid5 = test_shared_string5.getGroup();
    long int old_gid1 = test_shared_string1.getGroup();
    long int new_gid1 = old_gid1 + 10;
    gmx::replaceAllSharedStrings("TestString1 E", old_gid1, "TestString4 E", new_gid1);

    debug << "final strings after renaming via replaceAllSharedStrings(\"TestString1 E\"," << old_gid1 << ", \"TestString4 E\"," << new_gid1 << ")" << std::endl;
    debug << "test_shared_string1 = \"" << test_shared_string1 << "\", group = " << test_shared_string1.getGroup() << ", useCountP = " << test_shared_string1.getUseCountP() << ", useCountPP = " << test_shared_string1.getUseCountPP() << ", use_count = " << test_shared_string1.getUseCount() << std::endl;
    debug << "test_shared_string2 = \"" << test_shared_string2 << "\", group = " << test_shared_string2.getGroup() << ", useCountP = " << test_shared_string2.getUseCountP() << ", useCountPP = " << test_shared_string2.getUseCountPP() << ", use_count = " << test_shared_string2.getUseCount() << std::endl;
    debug << "test_shared_string3 = \"" << test_shared_string3 << "\", group = " << test_shared_string3.getGroup() << ", useCountP = " << test_shared_string3.getUseCountP() << ", useCountPP = " << test_shared_string3.getUseCountPP() << ", use_count = " << test_shared_string3.getUseCount() << std::endl;
    debug << "test_shared_string4 = \"" << test_shared_string4 << "\", group = " << test_shared_string4.getGroup() << ", useCountP = " << test_shared_string4.getUseCountP() << ", useCountPP = " << test_shared_string4.getUseCountPP() << ", use_count = " << test_shared_string4.getUseCount() << std::endl;
    debug << "test_shared_string5 = \"" << test_shared_string5 << "\", group = " << test_shared_string5.getGroup() << ", useCountP = " << test_shared_string5.getUseCountP() << ", useCountPP = " << test_shared_string5.getUseCountPP() << ", use_count = " << test_shared_string5.getUseCount() << std::endl;
    if (debug)
    {
        debug << "HERE :test_shared_string1" << std::endl;
        test_shared_string1.debugStats();
        debug << "HERE :test_shared_string2" << std::endl;
        test_shared_string2.debugStats();
        debug << "HERE :test_shared_string3" << std::endl;
        test_shared_string3.debugStats();
        debug << "HERE :test_shared_string4" << std::endl;
        test_shared_string4.debugStats();
        debug << "HERE :test_shared_string5" << std::endl;
        test_shared_string5.debugStats();
        debug << std::endl << "HERE symbol table topology:" << std::endl;
        test_shared_string5.debugSymtabStats();
    }

    EXPECT_TRUE(test_shared_string1.std_str() == "TestString4 E" &&
                test_shared_string2.std_str() == "TestString4 E" &&
                test_shared_string3.std_str() == "TestString4 E" &&
                test_shared_string4.std_str() == "TestString4 E" &&
                test_shared_string5.std_str() == "TestString1 E" &&
                test_shared_string1 == test_shared_string2 &&
                test_shared_string1 == test_shared_string3 &&
                test_shared_string1 == test_shared_string4)
    << "Test replaceAll() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect content." << std::endl;

    /*! useCount  : total number of shared_string instances referencing the same symbol / string content
        useCountP : total number of shared_ptr instances from the same set as
                    the shared_ptr of this shared_string referencing the same symbol,
                    can differ from useCount after renaming of a different string to
                    another string that also already existed in the symbol table
        useCountPP: total number of inner shared_ptrs referencing the same symbol
     */

    EXPECT_TRUE(test_shared_string1.getGroup() == new_gid1 && test_shared_string5.getGroup() == old_gid5)
    << "Test replaceAll() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect use group IDs. Symbol table appears corrupted." << std::endl;

    EXPECT_TRUE(test_shared_string1.getUseCount() == 4)
    << "Test replaceAll() with source and target \"string\"s as arguments failed. Shared_string(s) have incorrect use count. Symbol table appears corrupted." << std::endl;

}

TEST_F(SharedStringTest, ReadFromFileSingleThreaded)
{
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
            error << "SharedStringTest ReadFromFileSingleThreaded: failed to open input file, " << textfile.c_str() << std::endl;
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
            debug << "Number of unique string(s) read = " << shared_vec[shared_vec.size()-1].getNr() << std::endl;
            debug << "Size of shared_string infrastructure = " << sizeof(gmx::shared_string) << std::endl;
            debug << "Memory occupied by symbol table (estimated) = " << shared_vec[0].getTotalMem() << " byte" << std::endl;
        }

        if (shared_vec.size() != 0)
        {
            std::string tmpstr1 = "What_A_Terribly_Looooong_Name_For_A_Binding_Form";
            std::string tmpstr2("A_Shorter_Name");
            // remove preexisting instances of the strings to stick to the simple case of one instance only
            gmx::replaceAllSharedStrings(tmpstr1, "ABC");
            gmx::replaceAllSharedStrings(tmpstr2, "ABC");
            shared_vec.push_back(gmx::shared_string("What_A_Terribly_Looooong_Name_For_A_Binding_Form"));

            shared_vec[shared_vec.size()-1] = tmpstr1;
            ssize_t mem_before = shared_vec[shared_vec.size()-1].getTotalMem();
            gmx::replaceAllSharedStrings(tmpstr1, tmpstr2);
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

            ssize_t mem_after = shared_vec[shared_vec.size()-1].getTotalMem();
            EXPECT_TRUE( ((static_cast<ssize_t>(tmpstr1.length()) - static_cast<ssize_t>(tmpstr2.length())) * static_cast<ssize_t>(sizeof(char)))
                         >= (static_cast<ssize_t>(mem_after)        - static_cast<ssize_t>(mem_before)) )
            << "Memory book keeping upon renaming seems not to work. difference expected = "
            << ((static_cast<ssize_t>(tmpstr2.length()) - static_cast<ssize_t>(tmpstr1.length())) * static_cast<ssize_t>(sizeof(char)))
            << " got = " << (static_cast<ssize_t>(mem_after) - static_cast<ssize_t>(mem_before))
            << std::endl;

            // test symbol table cleanup
            // introduce two unique strings at the end, and remove the first one
            // to create a vacant position / gap in the symbol table for clean()
            shared_vec.push_back("Uniquetest_shared_stringXYZxyzABCabc-0");
            shared_vec.push_back("Uniquetest_shared_stringXYZxyzABCabc-1");
            shared_vec[shared_vec.size()-2] = shared_vec[0];
            mem_before = shared_vec[0].getTotalMem();
            shared_vec[0].clean();
            mem_after = shared_vec[0].getTotalMem();
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

}

//! this one is intended for testing with helgrind or similar for threading errors
//! with GNU libstdc++, setting the environment variable GLIBCPP_FORCE_NEW to 1
//! avoids false-positives in detecting data races due to reallocating memory in
//! regions previously allocated and freed again by other threads
TEST_F(SharedStringTest, ReadFromFileMultiThreaded)
{
    debug << "Test multi-threaded reading from text file into shared_strings  " << std::endl;
    debug << "and test some functions (insert/retrieve/remove/rename) entries," << std::endl;
    debug << "report (approximate) total amount of occupied memory) ..."        << std::endl;

    const size_t             nthreads = (std::thread::hardware_concurrency() > 1 ? std::thread::hardware_concurrency() : 2);
    std::vector<std::thread> tvec;
    debug << "Number of threads used = " << nthreads << std::endl;

    const size_t             ntextfiles = 3;
    debug << "Number of textfiles assigned round-robin to the reading threads = " << ntextfiles << std::endl;
    std::vector<std::string> textfiles;
    for (size_t i = 0; i < ntextfiles; ++i)
    {
        char buffer[50];
        sprintf(buffer, "refdata/textfile.%lu.txt", i + 1);
        textfiles.push_back(std::string(buffer));
    }
    // map ntextfiles input files to nthreads threads for tests with arbitrary numbers of threads,
    // regardless of the number of available text files
    bool                     all_resolved = true;
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
            debug << "Number of unique string(s) read = " << shared_vec[shared_vec.size()-1].getNr() << std::endl;
            debug << "Size of shared_string infrastructure = " << sizeof(gmx::shared_string) << std::endl;
            debug << "Memory occupied by symbol table (estimated) = " << shared_vec[0].getTotalMem() << " byte" << std::endl;
        }
    }

}


} // namespace
