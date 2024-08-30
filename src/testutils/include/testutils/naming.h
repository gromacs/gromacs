/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Functionality for customizing names of tests and the files
 * containing their reference data in a declarative way.
 *
 * To provide GoogleTest with helpfully customized names for test cases
 * constructed via testing::TestWithParam<TestParameters>,
 * - declare a functor of type
 *   \c gmx::test::NameOfTestFromTuple<TestParameters>,
 * - give its constructor formatting functions, functors, lambdas, or
 *   EnumerationArrays in a tuple whose elements correspond to
 *   the tuple that parameterizes the test, and
 * - pass that functor as the optional fourth parameter
 *   to INSTANTIATE_TEST_SUITE_P.
 * GoogleTest will use the functor and produce suitable names.
 * An example follows:
 *
 * \code
 *
 * //! Enumeration used as an element of the test parameter tuple
 * enum class : int Flavor { A, B, C, Count };
 * EnumerationArray<Flavor, const char*> sc_flavorName = { "A", "B", "C" };
 *
 * //! The test-parameter tuple
 * std::tuple<int, std::string, Flavor, float> TestParameters;
 * //! Parameterized test fixture
 * class ExampleTest : public ::testing::TestWithParam<TestParameters> { };
 *
 * //! Functor containing tuple for naming each test case
 * const NameOfTestFromTuple<TestParameters> sc_testNamer{
 *   std::make_tuple(PrefixFormatter<int, intToString>{"i_"},
 *                   useString,
 *                   sc_flavorName,
 *                   doubleToString) };
 *
 * // Extra functor for further customizing the filename used for the
 * // reference data, particularly useful when some parameters should
 * // not influence the name of the file. See below for details.
 * // This example makes reference-data filenames like
 * // Combinations_ExampleTest_i_0_foo_3_1416.xml
 * const RefDataFilenameMaker<TestParameters> sc_refDataFilenameMaker {
 *   std::make_tuple(PrefixFormatter<int, intToString>{"i_"},
 *                   useString,
 *                   toEmptyString, // ignores the Flavor
 *                   doubleToString) };
 *
 * TEST_P(ExampleTest, Works)
 * {
 *   // Use C++17 structured bindings to extract the test parameters
 *   auto [ intParam, stringParam, flavorParam, floatParam ] = GetParam();
 *   // Use the functor to name the reference data filename appropriately,
 *   // ie. where it is expected that Flavor does not change the values expected.
 *   TestReferenceData refData(sc_refDataFilenameMaker(GetParam()));
 *
 *   // Go test things against the reference data
 *   // ...
 * }
 *
 * // Use the functor to name the test cases
 * INSTANTIATE_TEST_SUITE_P(Combinations, ExampleTest,
 *                          ::testing::Combine(::testing::ValuesIn(0, 1, 3),
 *                                             ::testing::Values("foo", "bar"),
 *                                             ::testing::Values(Flavor::A, Flavor::B),
 *                                             ::testing::Values(3.1416, 2.7183)),
 *                          sc_testNamer);
 *
 * \endcode
 *
 * When passed, the \c sc_testNamer functor is called by GoogleTest to
 * generate a unique name for the test, which is done by \c
 * NameOfTestFromTuple<TestParamers> calling the formatters (passing
 * the respective parameters) to construct an underscore-separated
 * string as the test name.
 *
 * Several helper functions have been provided to smooth common
 * use cases, including
 * - \c useString for using a string paramater as-is
 * - \c PrefixFormatter for prefixing an existing formatter with a label string
 * - \c toEmptyString for producing no output for a string
 *
 * \c toEmptyString is particularly useful for the case where the test
 * case needs reference data, but only a subset of the test
 * parameters contribute to computing the expected result, while others
 * specify the implementation (e.g. target hardware, or which sets of
 * outputs are to be computed, or an implementation detail). To do this:
 * - declare an additional functor of type \c RefDataFilenameMaker<TestParameters>
 *   and provide formatters to it to either contribute to or be omitted from
 *   the name of the reference data file (mostly these are the same as
 *   the formatters used with \c sc_testNamer),
 * - construct the \c TestReferenceData with the name returned by
 *   `sc_refDataFilenameMaker(GetParam())`
 * - expect that the same reference data can be used transparently from
 *   multiple tests if some parameters do not contribute to the file naming.
 * The above example shows such a use case.
 *
 * Alternatively, if the test author prefers to use a struct for parameters,
 * the supporting machinery can be declared inside that struct for clarity,
 * like
 *
 * \code
 *
 * //! The test parameters
 * struct TestParameters
 * {
 *   //! The test-parameter tuple, must match the layout of the rest of the class
 *   using ParametersTuple = std::tuple<int, std::string, Flavor, float>;
 *   //! Some int
 *   int i;
 *   // ... other parameters
 *   //! Tuple of formatters to name the parameterized test cases
 *   static const NameOfTestFromTuple<ParametersTuple> sc_testNamer;
 *   //! Tuple of formatters to name the test-case reference data uniquely enough
 *   static const RefDataFilenameMaker<ParametersTuple> sc_refDataFilenameMaker;
 * };
 * // Then define the static variables outside the class in the same manner as
 * // above, as required by the language.
 *
 * \endcode
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_NAMING_H
#define GMX_TESTUTILS_NAMING_H

#include <algorithm>
#include <filesystem>
#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

#include <gtest/gtest.h>

#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/stringutil.h"

#include "testasserts.h"

namespace gmx
{

namespace test
{

namespace detail
{

/*! \libinternal \brief
 * Single-component variant for a formatting callable for variables of type \c T
 *
 * When the user wants to supply a formatting function for the name of
 * a GoogleTest parameters of an arbitrary non-enum type, a callable
 * taking a parameter of type \c T is required. This variant can
 * always contain it, whether it is a function, functor, or lambda. */
template<typename T, typename Enable = void>
struct FormatterVariant
{
    using Variant = std::variant<std::function<std::string(T)>>;
};

/*! \libinternal \brief
 * Specialization for formatting function of variables whose type is an enum
 *
 * When the user wants to use an enum as the type of a GoogleTest test
 * parameter, they must choose a class enum and can choose to supply
 * either a formatting callable or a matching EnumerationArray, and
 * this variant can contain either case. */
template<typename Enum>
struct FormatterVariant<Enum, typename std::enable_if_t<std::is_enum_v<Enum>>>
{
    static_assert(EnumClassSuitsEnumerationArray<Enum>::value,
                  "Enum parameter to GoogleTest test case must be suitable for EnumerationArray, "
                  "ie. have a Count field");
    using Variant = std::variant<std::function<std::string(Enum)>, EnumerationArray<Enum, const char*>>;
};

/*! \libinternal \brief
 * Template for mapping parameter types to formatter variants */
template<typename Tuple>
struct ParamsToFormatterVariants;

/*! \libinternal
 * \brief Specialization for tuple of parameter types to tuple of formatter variants
 *
 * This makes it easy to declare functions and classes that take
 * parameters that are tuples of parameters and matching tuples of
 * formatter variants. */
template<typename... Ts>
struct ParamsToFormatterVariants<std::tuple<Ts...>>
{
    using type = std::tuple<typename FormatterVariant<Ts>::Variant...>;
};

/*! \libinternal
 * \brief Type trait for whether \c T is a std::tuple
 *
 * There is no std::is_tuple<T>. */
template<typename T>
struct IsTuple : std::false_type
{
};

/*! \libinternal \brief
 * Specialization to successfully match a std::tuple */
template<typename... Ts>
struct IsTuple<std::tuple<Ts...>> : std::true_type
{
};

/*! \brief Helper constant for the std::visitor in \c formatNameFromParam
 *
 * This is a readable way to fail a static assertion while ensuring that
 * the name of a problematic type appears in the message. */
template<class>
inline constexpr bool sc_alwaysFalse = false;

/*! \brief Apply the \c formatterVariant to the respective \c param to produce a string
 *
 * The resulting string is intended to help name the part of a
 * GoogleTest test case that has the parameter value \c param.
 *
 * \tparam Param  The type of a parameter used to declare the
 *                test fixture class via testing::TestWithParam.
 */
template<typename Param>
std::string formatNameFromParam(const Param                                     param,
                                const typename FormatterVariant<Param>::Variant formatterVariant)
{
    // Use std::visit to use the \c formatterVariant on \c param to
    // obtain a string.
    return std::visit(
            [&](auto&& formatter) {
                using Formatter = std::decay_t<decltype(formatter)>;
                // If this is a function-like thing taking \c param, call it.
                if constexpr (std::is_invocable_v<Formatter, Param>)
                {
                    if (formatter == nullptr)
                    {
                        GMX_THROW(APIError("Must have valid formatter"));
                    }
                    return std::string{ formatter(param) };
                }
                // If formatter is an EnumerationArray indexed by an
                // enum param, look it up. But first do some
                // compile-time checks that provide useful compiler
                // error messages when something is mismatched.
                else if constexpr (!std::is_enum_v<Param>)
                {
                    static_assert(
                            sc_alwaysFalse<Formatter>,
                            "When formatter is not a callable, the parameter must be an enum");
                }
                else if constexpr (!EnumClassSuitsEnumerationArray<Param>::value)
                {
                    static_assert(sc_alwaysFalse<Formatter>,
                                  "Enum parameter to test case must be suitable for "
                                  "EnumerationArray (ie. have a Count field)");
                }
        // All nvcc until 11.8.0 and until at least 12.3.1 mis-compile
        // this if-constexpr nest by failing to short-circuit the
        // logic. This leads them to try to instantiate
        // EnumerationArray<Param, const char*> when Param is not an
        // enum type, which ought to have been excluded by the
        // conditions above. That instantiation fails, which causes a
        // compilation failure. Other compilers are fine with this
        // code, so we leave it as a developer convenience.
#if !defined(__CUDACC_VER_MAJOR__) || (__CUDACC_VER_MAJOR__ > 12) \
        || ((__CUDACC_VER_MAJOR__ == 12) && (__CUDACC_VER_MINOR__ > 3))
                else if constexpr (!std::is_same_v<Formatter, EnumerationArray<Param, const char*>>)
                {
                    static_assert(sc_alwaysFalse<Formatter>,
                                  "When formatter taking an enum is not a callable, it must be an "
                                  "EnumerationArray");
                }
#endif
                else
                {
                    // The type of \c param is an enumeration suitable
                    // for \c gmx::EnumerationArray and formatter is
                    // such a \c gmx::EnumerationArray.
                    return std::string{ formatter[param] };
                }
#if !defined(__CUDACC_VER_MAJOR__) || (__CUDACC_VER_MAJOR__ > 11) \
        || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ > 4))
        // Compiler recognizes that we can't get here and so
        // we don't need to write an unreachable return
        // statement.
#else
                // nvcc before 11.5.0 fails to recognize that the above
                // if-constexpr nest always produced a return statement
                // so complains unless we add this unnecessary one.
                return std::string{};
#endif
            },
            formatterVariant);
}

/*! \brief Fold a constexpr sequence of integers that match \c S to calls
 * to the function \c f taking a matching std::integral_constant. */
template<typename T, T... S, typename F>
constexpr void ForSequence(std::integer_sequence<T, S...>, F&& f)
{
    (void(f(std::integral_constant<T, S>{})), ...);
}

/*! \brief Apply the \c formatters to the respective \c params to
 * produce an underscore-separated string
 *
 * The resulting string is intended to help name a GoogleTest test
 * case.
 *
 * If a formatter produces an empty string, no underscore will be
 * emitted following it.
 *
 * \tparam ParametersTuple  A tuple of test parameters used to declare the
 *                          test fixture class via testing::TestWithParam.
 */
template<typename ParametersTuple>
std::string mapNameFormattersToParameters(const ParametersTuple params,
                                          const typename ParamsToFormatterVariants<ParametersTuple>::type& formatters)
{
    std::vector<std::string> testNameComponents;
    constexpr size_t         sizeOfParametersTuple = std::tuple_size_v<ParametersTuple>;
    constexpr size_t         sizeOfFormatterVariantsTuple =
            std::tuple_size_v<typename ParamsToFormatterVariants<ParametersTuple>::type>;
    static_assert(sizeOfParametersTuple == sizeOfFormatterVariantsTuple,
                  "Cannot combine parameter and formatter tuples of different sizes");
    // Do a compile-time "iteration" over the elements of \c param and
    // \c formatters to apply an element of the latter to the
    // corresponding element of the former to produce a chunk of the
    // test name.
    ForSequence(std::make_index_sequence<sizeOfParametersTuple>{}, [&](auto i) {
        // Here we are in fact using the local variable `i` as
        // constexpr, even though it can't be declared that way
        // because it is actually a std::integral_constant<size_t, N>
        // for some N from the sequence.
        const auto param                        = std::get<i>(params);
        using Param                             = std::decay_t<decltype(param)>;
        using FormatterVariant                  = typename FormatterVariant<Param>::Variant;
        const FormatterVariant formatterVariant = std::get<i>(formatters);
        // Call a helper function to actually apply \c
        // formatterVariant to \c param.
        const std::string result = formatNameFromParam(param, formatterVariant);
        if (!result.empty())
        {
            testNameComponents.push_back(result);
        }
    });
    std::string testName = joinStrings(testNameComponents, "_");
    // Note that the returned name must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");
    testName = replaceAll(testName, "/", "_");
    return testName;
}

} // namespace detail

/*! \libinternal
 * \brief Function object that helps GoogleTest name our test cases
 *
 * \tparam ParametersTuple  A tuple of test parameters used to declare the
 *                          test fixture class via testing::TestWithParam.
 */
template<typename ParametersTuple>
class NameOfTestFromTuple
{
public:
    static_assert(detail::IsTuple<ParametersTuple>::value);
    using Formatters = typename detail::ParamsToFormatterVariants<ParametersTuple>::type;

    //! Constructor
    NameOfTestFromTuple(Formatters formatters) : formatters_(formatters) {}
    /*! \brief Return the name of this test.
     *
     * Called by GoogleTest while instantiating value-parameterized
     * tests with the test \c info that contains a tuple of parameters
     * which must be converted to a unique name. */
    std::string operator()(const testing::TestParamInfo<ParametersTuple>& info) const
    {
        using InputParameters = decltype(info.param);
        static_assert(detail::IsTuple<InputParameters>::value);
        return detail::mapNameFormattersToParameters(info.param, formatters_);
    }

private:
    /*! \brief A tuple of \c FormatterVariant objects that match \c
     * ParametersTuple.
     *
     * The formatters convert the respective values in a \c
     * ParametersTuple to unique strings to help build a unique
     * name for the test case. */
    Formatters formatters_;
};

/*! \libinternal \brief
 * Functor to name a refdata filename for this test
 *
 * We sometimes want the same reference data to apply to multiple test
 * cases, e.g. because we implement multiple ways to do the same thing
 * and want to check that all implementations match the same reference
 * data.
 *
 * That means we need to store the reference data in a file whose name
 * relates only to the test parameters that change the values, and not
 * those that describe different implementations. In such cases, we
 * cannot simple re-use the name provided by \c NameOfTestFromTuple.
 *
 * By default, the reference data filename is set via a call to
 * gmx::TestFileManager::getTestSpecificFileName() that queries
 * GoogleTest and gets a string that includes the return value for
 * NameOfTestFromTuple::operator(). This code works similarly, but
 * generally the user passes different formatters for some parameters,
 * so that some of the test parameters do not contribute to the name
 * produced. */
template<typename ParametersTuple>
class RefDataFilenameMaker
{
public:
    using Formatters = typename detail::ParamsToFormatterVariants<ParametersTuple>::type;

    //! Constructor
    RefDataFilenameMaker(Formatters formatters) : formatters_(formatters) {}
    //! Functor to make the reference-data filename requested
    std::filesystem::path operator()(const ParametersTuple& params) const
    {
        static_assert(detail::IsTuple<ParametersTuple>::value);
        std::string testParametersName = detail::mapNameFormattersToParameters(params, formatters_);

        // Build the complete filename like
        // gmx::TestFileManager::getTestSpecificFilename() does it, but
        // there's no requirement that they match.
        const ::testing::TestInfo* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string                testSuiteName(testInfo->test_suite_name());
        std::string                testName(testInfo->name());
        // If the test fixture body was created with TEST_P then there is
        // a further user-provided name that should contribute to the
        // reference data name. If present, it precedes "/".
        auto        separatorPos = testName.find("/");
        std::string testFixtureBodyName =
                (separatorPos == std::string::npos) ? "" : (testName.substr(0, separatorPos));
        std::string refDataName = testSuiteName;
        // Add separators only where needed
        if (!testFixtureBodyName.empty())
        {
            refDataName += "_";
        }
        refDataName += testFixtureBodyName;
        if (!testParametersName.empty())
        {
            refDataName += "_";
        }
        refDataName += testParametersName;
        // Filenames shouldn't have '/' characters in them (for
        // sanity).
        std::replace(refDataName.begin(), refDataName.end(), '/', '_');
        // Check that the name isn't too long
        checkTestNameLength(refDataName);
        // Turn the stem into a filename.
        return refDataName + ".xml";
    }

private:
    /*! \brief A tuple of \c FormatterVariant objects that match \c
     * ParametersTuple.
     *
     * The formatters convert the respective values in a \c
     * ParametersTuple to unique strings to help build a unique
     * enough name for the reference data for the test case. */
    const Formatters formatters_;
};

//! Formatter to pass std::string through
static inline std::string useString(const std::string s)
{
    return s;
}

/*! \libinternal
 * \brief Functor to add a prefix to the return from \c formatter for \c T
 */
template<typename T, std::string (*formatter)(T)>
struct PrefixFormatter
{
    std::string prefix;
    std::string operator()(const T& t) { return prefix + formatter(t); }
};

//! Formatter to ignore test parameters in a self-documenting way
template<typename T>
std::string inline toEmptyString(const T /* t */)
{
    return std::string{};
}

} // namespace test
} // namespace gmx

#endif
