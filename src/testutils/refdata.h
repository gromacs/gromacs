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
/*! \libinternal \file
 * \brief
 * Functionality for writing tests that can produce their own reference data.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_REFDATA_H
#define GMX_TESTUTILS_REFDATA_H

#include <string>

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Mode of operation for reference data handling.
 *
 * There should be no need to use this type outside the test utility module.
 */
enum ReferenceDataMode
{
    /*! \brief
     * Compare to existing reference data.
     *
     * If reference data does not exist, or if the test results differ from
     * those in the reference data, the test fails.
     */
    erefdataCompare,
    /*! \brief
     * Create missing reference data.
     *
     * If reference data does not exist for a test, that test behaves as if
     * ::erefdataUpdateAll had been specified.  Tests for which reference data
     * exists, behave like with ::erefdataCompare.
     */
    erefdataCreateMissing,
    /*! \brief
     * Update reference data, overwriting old data.
     *
     * Tests utilizing reference data should always pass in this mode unless
     * there is an I/O error.
     */
    erefdataUpdateAll
};

/*! \libinternal \brief
 * Returns the global reference data mode.
 *
 * There should be no need to use this function outside the test utility module.
 */
ReferenceDataMode getReferenceDataMode();
/*! \libinternal \brief
 * Sets the global reference data mode.
 *
 * There should be no need to use this function outside the test utility module.
 */
void setReferenceDataMode(ReferenceDataMode mode);
/*! \libinternal \brief
 * Returns the directory where reference data files are stored.
 *
 * There should be no need to use this function outside the test utility module.
 */
std::string getReferenceDataPath();
/*! \libinternal \brief
 * Initializes reference data handling.
 *
 * Sets the reference data mode based on command-line arguments.  By default,
 * ::erefdataCompare is used, but \c --create-ref-data or \c --update-ref-data
 * can be used to change it.
 * Recognized command-line arguments are removed from the list.
 *
 * Also creates the directory for storing the reference data.
 *
 * This function is automatically called by test_main_gtest.cpp and
 * test_main_gmock.cpp.
 */
int initReferenceData(int *argc, char **argv);

class TestReferenceChecker;

/*! \libinternal \brief
 * Handles creation of and comparison to test reference data.
 *
 * This class provides functionality to use the same code to generate reference
 * data and then on later runs compare the results of the code against that
 * reference.  The mode in which the class operates (writing reference data or
 * comparing against existing data) is set with parseReferenceDataArgs(), which
 * is automatically called when using the testutils module to implement tests.
 * Tests only need to create an instance of TestReferenceData, obtain a
 * TestReferenceChecker using the rootChecker() method and use the various
 * check*() methods in TestReferenceChecker to indicate values to check.  If
 * the test is running in reference data creation mode, it will produce an XML
 * file with the values recorder.  In comparison mode, it will read that same
 * XML file and produce a Google Test non-fatal assertion for every discrepancy
 * it detects with the reference data (including missing reference data file or
 * individual item).  Exceptions derived from TestException are thrown for I/O
 * errors and syntax errors in the reference data.
 *
 * Simple example (using Google Test):
 * \code
int functionToTest(int param);

TEST(MyTest, SimpleTest)
{
    gmx::test::TestReferenceData data;

    gmx::test::TestReferenceChecker checker(data.rootChecker());
    checker.checkInteger(functionToTest(3), "ValueWith3");
    checker.checkInteger(functionToTest(5), "ValueWith5");
    gmx::test::TestReferenceChecker compound(checker.startCompound("CustomCompound", "Item"));
    compound.checkInteger(function2ToTest(3), "ValueWith3");
    compound.checkInteger(function2ToTest(5), "ValueWith5");
    checker.checkInteger(functionToTest(4), "ValueWith4");
}
 * \endcode
 *
 * This class is only available if both Google Test and libxml2 are enabled.
 * If either one is missing, trying to use this class will result in unresolved
 * symbols in linking.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class TestReferenceData
{
    public:
        /*! \brief
         * Initializes the reference data in the global mode.
         */
        TestReferenceData();
        /*! \brief
         * Initializes the reference data in a specific mode.
         *
         * This function is mainly useful for self-testing the reference data
         * framework.
         * The default constructor should be used in tests utilizing this class.
         */
        explicit TestReferenceData(ReferenceDataMode mode);
        /*! \brief
         * Frees reference data structures.
         *
         * In the current implementation, this function writes the reference
         * data out if necessary.
         */
        ~TestReferenceData();

        //! Returns true if reference data is currently being written.
        bool isWriteMode() const;

        /*! \brief
         * Returns a root-level checker object for comparisons.
         *
         * Each call returns an independent instance.
         */
        TestReferenceChecker rootChecker();

    private:
        class Impl;

        Impl                   *_impl;

        // Disallow copy and assign.
        TestReferenceData(const TestReferenceData &);
        void operator =(const TestReferenceData &);
};

/*! \libinternal \brief
 * Handles comparison to test reference data.
 *
 * Every check*() method takes an id string as the last parameter.  This id is
 * used to uniquely identify the value in the reference data, and it makes the
 * output XML more human-friendly and more robust to errors.  The id can be
 * NULL; in this case, multiple elements with no id are created, and they will
 * be matched in the same order as in which they are created.  The
 * checkCompound() method can be used to create a set of reference values
 * grouped together.  In this case, all check*() calls using the returned child
 * TestReferenceChecker object will create the reference data within this
 * group, and the ids only need to be unique within the compound.  Compounds
 * can be nested.
 *
 * For usage example, see TestReferenceData.
 *
 * Copies of this class behave have independent internal state.
 *
 * This class is only available if both Google Test and libxml2 are enabled.
 * If either one is missing, trying to use this class will result in unresolved
 * symbols in linking.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class TestReferenceChecker
{
    public:
        /*! \brief
         * Creates a deep copy of the other checker.
         */
        TestReferenceChecker(const TestReferenceChecker &other);
        ~TestReferenceChecker();

        TestReferenceChecker &operator =(const TestReferenceChecker &other);

        //! Returns true if reference data is currently being written.
        bool isWriteMode() const;

        /*! \brief
         * Initializes comparison of a group of related data items.
         *
         * \param[in] type Informational type for the compound.
         * \param[in] id   Unique identifier for the compound among its
         *                 siblings.
         * \returns   Checker to use for comparison within the compound.
         *
         * All checks performed with the returned checker only
         * need to have unique ids within the compound, not globally.
         *
         * Compound structures can be nested.
         */
        TestReferenceChecker checkCompound(const char *type, const char *id);

        //! Check a single boolean value.
        void checkBoolean(bool value, const char *id);
        //! Check a single string value.
        void checkString(const char *value, const char *id);
        //! Check a single string value.
        void checkString(const std::string &value, const char *id);
        //! Check a single integer value.
        void checkInteger(int value, const char *id);
        //! Check a single single-precision floating point value.
        void checkFloat(float value, const char *id);
        //! Check a single double-precision floating point value.
        void checkDouble(double value, const char *id);
        //! Check a single floating point value.
        void checkReal(float value, const char *id);
        //! Check a single floating point value.
        void checkReal(double value, const char *id);
        //! Check a vector of three integer values.
        void checkVector(int value[3], const char *id);
        //! Check a vector of three single-precision floating point values.
        void checkVector(float value[3], const char *id);
        //! Check a vector of three double-precision floating point values.
        void checkVector(double value[3], const char *id);
        //! Check an arbitrary-length sequence of integer values.
        void checkSequenceInteger(size_t length, int *values,
                                  const char *id);
        //! Check an arbitrary-length sequence of double-precision values.
        void checkSequenceDouble(size_t length, double *values,
                                 const char *id);
        //! Check an arbitrary-length sequence of integer vector values.
        void checkSequenceVector(size_t length, int values[][3],
                                 const char *id);
        //! Check an arbitrary-length sequence of floating-point vector values.
        void checkSequenceVector(size_t length, float values[][3],
                                 const char *id);
        //! Check an arbitrary-length sequence of floating-point vector values.
        void checkSequenceVector(size_t length, double values[][3],
                                 const char *id);

    private:
        class Impl;

        /*! \brief
         * Constructs a checker with a specific internal state.
         *
         * Is private to only allow users of this class to create instances
         * using TestReferenceData::rootChecker() or checkCompound()
         * (or by copying).
         */
        explicit TestReferenceChecker(Impl *impl);

        Impl                   *_impl;

        /*! \brief
         * Needed to expose the constructor only to TestReferenceData.
         */
        friend class TestReferenceData;
};

} // namespace test
} // namespace gmx

#endif
