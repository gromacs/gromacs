/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Functionality for writing tests that can produce their own reference data.
 *
 * See \ref page_refdata for more details.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_REFDATA_H
#define GMX_TESTUTILS_REFDATA_H

#include <iterator>
#include <memory>
#include <string>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class IOptionsContainer;
class KeyValueTreeObject;
class KeyValueTreeValue;
class Any;

//! \internal
template<typename ValueType>
class BasicVector;

namespace test
{

class FloatingPointTolerance;

/*! \libinternal \brief
 * Mode of operation for reference data handling.
 *
 * There should be no need to use this type outside the test utility module.
 *
 * \ingroup module_testutils
 */
enum class ReferenceDataMode : int
{
    /*! \brief
     * Compare to existing reference data.
     *
     * If reference data does not exist, or if the test results differ from
     * those in the reference data, the test fails.
     */
    Compare,
    /*! \brief
     * Create missing reference data.
     *
     * If reference data does not exist for a test, that test behaves as if
     * ReferenceDataMode::UpdateAll had been specified.  Tests for which reference data
     * exists, behave like with ReferenceDataMode::Compare.
     */
    CreateMissing,
    /*! \brief
     * Update reference data that does not pass comparison.
     *
     * Tests utilizing reference data should always pass in this mode unless
     * there is an I/O error.
     */
    UpdateChanged,
    /*! \brief
     * Update reference data, overwriting old data.
     *
     * Tests utilizing reference data should always pass in this mode unless
     * there is an I/O error.
     */
    UpdateAll,
    //! Marks the end of the enum
    Count
};

/*! \libinternal \brief
 * Returns the global mode describing how reference data is handled. */
ReferenceDataMode referenceDataMode();

/*! \libinternal \brief
 * Initializes reference data handling.
 *
 * Adds command-line options to \p options to set the reference data mode.
 * By default, ReferenceDataMode::Compare is used, but ``--ref-data create`` or
 * ``--ref-data update`` can be used to change it.
 *
 * This function is automatically called by initTestUtils().
 *
 * \ingroup module_testutils
 */
void initReferenceData(IOptionsContainer* options);

class TestReferenceChecker;

namespace internal
{
class TestReferenceDataImpl;
}

/*! \libinternal \brief
 * Handles creation of and comparison to test reference data.
 *
 * See \ref page_refdata for an overview of the functionality.
 *
 * This class provides functionality to use the same code to generate reference
 * data and then on later runs compare the results of the code against that
 * reference.  The mode in which the class operates (writing reference data or
 * comparing against existing data) is set using a command-line option that
 * is automatically managed when using the testutils module to implement tests.
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

   namespace gmx
   {
   namespace test
   {
   TEST(MyTest, SimpleTest)
   {
       TestReferenceData data;

       TestReferenceChecker checker(data.rootChecker());
       checker.checkInteger(functionToTest(3), "ValueWith3");
       checker.checkInteger(functionToTest(5), "ValueWith5");
       TestReferenceChecker compound(
               checker.checkCompound("CustomCompound", "Item"));
       compound.checkInteger(function2ToTest(3), "ValueWith3");
       compound.checkInteger(function2ToTest(5), "ValueWith5");
       checker.checkInteger(functionToTest(4), "ValueWith4");
       checker.checkVector(functionProducingRVec(), "Describe The RVec");
   }
   } // namespace test
   } // namespace gmx
 * \endcode
 *
 * If rootChecker() is never called, no comparison is done (i.e., missing
 * reference data file is not reported as an error, nor is empty reference data
 * file created in write mode).
 *
 * For floating-point comparisons, the reference data should be generated in
 * double precision (currently, no warning is provided even if this is not the
 * case, but the double precision tests will then very likely fail).
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
    /*! \brief Initializes the reference data in the global mode while
     * providing the test name.
     *
     * Providing a name for the test, rather than using the default
     * name generated by GoogleTest, is useful when many tests should
     * share common reference data because e.g. they are
     * value-parameterized over multiple equivalent implementations or
     * run modes. The name should include a filename suffix
     * e.g. `.xml` and refers to a file that will be found in the
     * `refdata` subdirectory for the test case. */
    explicit TestReferenceData(std::string testNameOverride);
    /*! \brief
     * Initializes the reference data in a specific mode.
     *
     * This function is only useful for self-testing the reference data
     * framework.  As such, it also puts the framework in a state where it
     * logs additional internal information for failures to help diagnosing
     * problems in the framework, and stores the reference data in a
     * temporary directory instead of the source tree.
     * The default constructor should be used in tests utilizing this class.
     */
    explicit TestReferenceData(ReferenceDataMode mode);
    /*! \brief
     * Frees reference data structures.
     *
     * The reference data is written out if necessary automatically when
     * the test finishes.
     */
    ~TestReferenceData();

    /*! \brief
     * Returns a root-level checker object for comparisons.
     *
     * Each call returns an independent instance.
     */
    TestReferenceChecker rootChecker();

private:
    std::shared_ptr<internal::TestReferenceDataImpl> impl_;

    GMX_DISALLOW_COPY_AND_ASSIGN(TestReferenceData);
};

/*! \libinternal \brief
 * Handles comparison to test reference data.
 *
 * Every check*() method takes an id string as th
 * e last parameter.  This id is
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
 * \inlibraryapi
 * \ingroup module_testutils
 */
class TestReferenceChecker
{
public:
    /*! \brief
     * Creates a checker that cannot be used for checking.
     *
     * Attempting to call the check methods generates an assert.
     * It is possible to check whether the checker is initialized by
     * calling isValid().
     * This constructor exists to allow declaring checker variables that
     * will receive their value later without resorting to dynamic
     * allocation.
     */
    TestReferenceChecker();
    //! Creates a deep copy of the other checker.
    explicit TestReferenceChecker(const TestReferenceChecker& other);
    //! Moves the checker.
    TestReferenceChecker(TestReferenceChecker&& other) noexcept;
    ~TestReferenceChecker();

    //! Prevents implicit copying during assignment.
    TestReferenceChecker& operator=(const TestReferenceChecker&) = delete;
    //! Assigns a test reference checker.
    TestReferenceChecker& operator=(TestReferenceChecker&& other) noexcept;

    //! Returns whether the checker is initialized.
    bool isValid() const;
    //! Allows testing whether the checker is initialized directly with if.
    explicit operator bool() const { return isValid(); }

    /*! \brief
     * Sets the tolerance for floating-point comparisons.
     *
     * All following floating-point comparisons using this checker will use
     * the new tolerance.  Child checkers created with checkCompound()
     * will inherit the tolerance from their parent checker at the time
     * checkCompound() is called.
     *
     * Does not throw.
     */
    void setDefaultTolerance(const FloatingPointTolerance& tolerance);

    /*! \brief
     * Checks that all reference values have been compared against.
     *
     * All values under the compound represented by this checker are
     * checked, and a non-fatal Google Test assertion is produced if some
     * values have not been used.
     *
     * If not called explicitly, the same check will be done for all
     * reference data values when the test ends.
     *
     * This method also marks the values used, so that subsequent checks
     * (including the check at the end of the test) will not produce
     * another assertion about the same values.
     */
    void checkUnusedEntries();

    /*! \brief Disables checking for unused entries
     *
     * \see checkUnusedEntries()
     */
    void disableUnusedEntriesCheck();

    /*! \brief
     * Checks whether a data item is present.
     *
     * \param[in] bPresent  Whether to check for presence or absence.
     * \param[in] id        Unique identifier of the item to check.
     * \returns   true if bPresent was true and the data item was found.
     *
     * If \p bPresent is true, checks that a data item with \p id is
     * present, otherwise checks that the data item is absent.
     * If the check fails, a non-fatal Google Test assertion is generated.
     *
     * If reference data is being written, the check always succeeds and the
     * return value is \p bPresent.
     *
     * The main use of this method is to assign meaning for missing
     * reference data.  Example use:
     * \code
       if (checker.checkPresent(bHaveVelocities, "Velocities"))
       {
           // <check the velocities>
       }
     * \endcode
     */
    bool checkPresent(bool bPresent, const char* id);

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
    TestReferenceChecker checkCompound(const char* type, const char* id);
    //! \copydoc checkCompound(const char *, const char *)
    TestReferenceChecker checkCompound(const char* type, const std::string& id);

    //! Check a single boolean value.
    void checkBoolean(bool value, const char* id);
    //! Check a single string value.
    void checkString(const char* value, const char* id);
    //! Check a single string value.
    void checkString(const std::string& value, const char* id);
    /*! \brief
     * Check a multi-line string value.
     *
     * This method works as checkString(), but should be used for long
     * strings that may contain, e.g., newlines.  Typically used to check
     * formatted output, and attempts to make the output XML such that it
     * is easier to edit by hand to set the desired output formatting.
     */
    void checkTextBlock(const std::string& value, const char* id);
    //! Check a single char value.
    void checkChar(char value, const char* id);
    //! Check a single unsigned char value.
    void checkUChar(unsigned char value, const char* id);
    //! Check a single integer value.
    void checkInteger(int value, const char* id);
    //! Check a single int32 value.
    void checkInt32(int32_t value, const char* id);
    //! Check a single uint32 value.
    void checkUInt32(uint32_t value, const char* id);
    //! Check a single int64 value.
    void checkInt64(int64_t value, const char* id);
    //! Check a single uint64 value.
    void checkUInt64(uint64_t value, const char* id);
    //! Check a single single-precision floating point value.
    void checkFloat(float value, const char* id);
    //! Check a single double-precision floating point value.
    void checkDouble(double value, const char* id);
    //! Check a single floating point value.
    void checkReal(float value, const char* id);
    //! Check a single floating point value.
    void checkReal(double value, const char* id);
    //! Check a vector of three integer values.
    void checkVector(const int value[3], const char* id);
    //! Check a vector of three single-precision floating point values.
    void checkVector(const float value[3], const char* id);
    //! Check a vector of three double-precision floating point values.
    void checkVector(const double value[3], const char* id);
    //! Check a BasicVector of ints, ie. IVec
    void checkVector(const BasicVector<int>& value, const char* id);
    //! Check a BasicVector of floats, ie. RVec
    void checkVector(const BasicVector<float>& value, const char* id);
    //! Check a BasicVector of doubles, ie. DVec
    void checkVector(const BasicVector<double>& value, const char* id);
    //! Check a single floating-point value from a string.
    void checkRealFromString(const std::string& value, const char* id);
    //! Checks a any value that contains a supported simple type.
    void checkAny(const Any& value, const char* id);
    //! Checks a key-value tree rooted at a object.
    void checkKeyValueTreeObject(const KeyValueTreeObject& tree, const char* id);
    //! Checks a generic key-value tree value.
    void checkKeyValueTreeValue(const KeyValueTreeValue& value, const char* id);

    /*! \name Methods to read values from reference data
     *
     * These methods assume that a value with the given `id` has already
     * been created in the test with `check*()` methods, and that it has
     * the correct type.
     *
     * Currently, these methods do not work correctly if the reference data
     * file does not exist, so a test using them may fail with exceptions
     * before the reference data has been generated.
     * \{
     */
    //! Reads an unsigned char value.
    unsigned char readUChar(const char* id);
    //! Reads an integer value.
    int readInteger(const char* id);
    //! Reads a 32-bit integer value.
    int32_t readInt32(const char* id);
    //! Reads a 64-bit integer value.
    int64_t readInt64(const char* id);
    //! Reads a float value.
    float readFloat(const char* id);
    //! Reads a double value.
    double readDouble(const char* id);
    //! Reads a string value.
    std::string readString(const char* id);
    //! \}

    /*! \name Overloaded versions of simple checker methods
     *
     * These methods provide overloads under a single name for all the
     * methods checkBoolean(), checkString(), checkReal() and checkVector().
     * They are provided mainly to allow template implementations (such as
     * checkSequence()).  Typically callers should use the individually
     * named versions for greater clarity.
     * \{
     */
    //! Check a single boolean value.
    void checkValue(bool value, const char* id) { checkBoolean(value, id); }
    //! Check a single string value.
    void checkValue(const char* value, const char* id) { checkString(value, id); }
    //! Check a single string value.
    void checkValue(const std::string& value, const char* id) { checkString(value, id); }
    //! Check a single signed integer value
    void checkValue(int value, const char* id) { checkInteger(value, id); }
    //! Check a single signed integer value of width 64 bits.
    void checkValue(int64_t value, const char* id) { checkInt64(value, id); }
    //! Check a single unsigned integer value of width 64 bits.
    void checkValue(uint64_t value, const char* id) { checkUInt64(value, id); }
    //! Check a single single-precision floating point value.
    void checkValue(float value, const char* id) { checkFloat(value, id); }
    //! Check a single double-precision floating point value.
    void checkValue(double value, const char* id) { checkDouble(value, id); }
    //! Check a vector of three integer values.
    void checkValue(const int value[3], const char* id) { checkVector(value, id); }
    //! Check a vector of three single-precision floating point values.
    void checkValue(const float value[3], const char* id) { checkVector(value, id); }
    //! Check a vector of three double-precision floating point values.
    void checkValue(const double value[3], const char* id) { checkVector(value, id); }
    //! Check a BasicVector of integer values, ie. IVec.
    void checkValue(const BasicVector<int>& value, const char* id) { checkVector(value, id); }
    //! Check a BasicVector of float values, ie. RVec.
    void checkValue(const BasicVector<float>& value, const char* id) { checkVector(value, id); }
    //! Check a BasicVector of double values, ie. DVec.
    void checkValue(const BasicVector<double>& value, const char* id) { checkVector(value, id); }
    //! Check a generic key-value tree value.
    void checkValue(const KeyValueTreeValue& value, const char* id)
    {
        checkKeyValueTreeValue(value, id);
    }
    //! Check a single complex number value.
    void checkValue(const t_complex& value, const char* id)
    {
        checkValue(value.re, id);
        checkValue(value.im, id);
    }

    /*!\}*/

    /*! \brief
     * Generic method to check a sequence of simple values.
     *
     * \tparam Iterator  Input iterator that allows multiple (two) passes.
     *      Value type must be one of those accepted by checkValue(), or
     *      implicitly convertible to one.
     * \param[in] begin  Iterator to the start of the range to check.
     * \param[in] end    Iterator to the end of the range to check.
     * \param[in] id     Unique identifier for the sequence among its
     *                   siblings.
     */
    template<class Iterator>
    void checkSequence(Iterator begin, Iterator end, const char* id)
    {
        typename std::iterator_traits<Iterator>::difference_type length = std::distance(begin, end);
        TestReferenceChecker compound(checkSequenceCompound(id, length));
        for (Iterator i = begin; i != end; ++i)
        {
            compound.checkValue(*i, nullptr);
        }
    }
    /*! \brief
     * Generic method to check a sequence of custom values.
     *
     * \tparam Iterator    Input iterator that allows multiple (two) passes.
     * \tparam ItemChecker Functor to check an individual value. Signature
     *      void(TestReferenceChecker *, const T &), where T is the value
     *      type of \p Iterator.
     * \param[in] begin  Iterator to the start of the range to check.
     * \param[in] end    Iterator to the end of the range to check.
     * \param[in] id     Unique identifier for the sequence among its
     *                   siblings.
     * \param[in] checkItem  Functor to check an individual item.
     *
     * This method creates a compound checker \c compound within which all
     * values of the sequence are checked.  Calls \c checkItem(&compound, *i)
     * with that compound for each iterator \c i in the range [begin, end).
     * \p checkItem should use the various check methods in the passed
     * checker to check each value.
     *
     * This method can be used to check a sequence made of compound types.
     * Typically \p checkItem will create a compound within the passed
     * checker to check different aspects of the value that was passed
     * to it. Either NULL or a unique identifier string must be used for
     * the id value of that compound. */
    template<class Iterator, class ItemChecker>
    void checkSequence(Iterator begin, Iterator end, const char* id, ItemChecker checkItem)
    {
        typename std::iterator_traits<Iterator>::difference_type length = std::distance(begin, end);
        TestReferenceChecker compound(checkSequenceCompound(id, length));
        for (Iterator i = begin; i != end; ++i)
        {
            checkItem(&compound, *i);
        }
    }
    /*! \brief
     * Check an array of values.
     *
     * \tparam T  Type of values to check. Should be one of those accepted
     *      by checkValue(), or implicitly convertible to one.
     *
     * \param[in] length  Number of values to check.
     * \param[in] values  Pointer to the first value to check.
     * \param[in] id     Unique identifier for the sequence among its
     *                   siblings.
     *
     * This is a convenience method that delegates all work to
     * checkSequence().
     */
    template<typename T>
    void checkSequenceArray(size_t length, const T* values, const char* id)
    {
        checkSequence(values, values + length, id);
    }
    /*! \brief
     * Convenience method for checking that a sequence is empty.
     *
     * \param[in] id     Unique identifier for the sequence among its
     *                   siblings.
     *
     * This method provides a convenient solution for a case where there is
     * implicitly a sequence to be checked, but there is no pointer
     * available to the values since the sequence is empty.
     * Since this method does not require the type of the values, it can be
     * used in such cases easily.
     */
    void checkEmptySequence(const char* id);
    /*! \brief
     * Initializes a compound for a sequence of items.
     *
     * \param[in] id     Unique identifier for the sequence among its
     *                   siblings.
     * \param[in] length Number of items that will be in the sequence.
     * \returns   Checker to use for comparison within the sequence.
     *
     * This method can be used to check custom sequences where
     * checkSequence() is not appropriate.
     */
    TestReferenceChecker checkSequenceCompound(const char* id, size_t length);

private:
    class Impl;

    /*! \brief
     * Constructs a checker with a specific internal state.
     *
     * Is private to only allow users of this class to create instances
     * using TestReferenceData::rootChecker() or checkCompound()
     * (or by copying).
     */
    explicit TestReferenceChecker(Impl* impl);

    std::unique_ptr<Impl> impl_;

    /*! \brief
     * Needed to expose the constructor only to TestReferenceData.
     */
    friend class TestReferenceData;
};

} // namespace test
} // namespace gmx

#endif
