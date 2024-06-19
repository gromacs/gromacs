/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/mdtypes/checkpointdata.h"

#include <climits>
#include <cstdint>

#include <algorithm>
#include <array>
#include <filesystem>
#include <functional>
#include <iterator>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/real.h"

#include "testutils/testfilemanager.h"

namespace gmx::test
{
namespace
{

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Struct allowing to check if type is vector of serializable data
 */
//! \{
template<class T>
struct IsVectorOfSerializableType
{
    static bool const value = false;
};
template<class T>
struct IsVectorOfSerializableType<std::vector<T>>
{
    static bool const value = IsSerializableType<T>::value;
};
//! \}

/*! \internal
 * \brief Unified looping over test data
 *
 * This class allows to write a loop over test data as
 *     for (const auto& value : TestValues::testValueGenerator<type>())
 * where type can be any of std::string, int, int64_t, bool, float, double,
 * std::vector<[std::string, int, int64_6, float, double]>, or tensor.
 */
class TestValues
{
public:
    /*! \internal
     * \brief Helper class allowing to loop over test values
     * \tparam T  type of value
     */
    template<typename T>
    class TestValueGenerator
    {
    public:
        //! Custom iterator
        class Iterator
        {
        public:
            explicit Iterator(const T* ptr) : ptr_(ptr) {}
            Iterator operator++();
            bool     operator!=(const Iterator& other) const { return ptr_ != other.ptr_; }
            const T& operator*() const { return *ptr_; }

        private:
            const T* ptr_;
        };

        Iterator begin() const;
        Iterator end() const;
    };

    /*! \internal
     * \brief Static function returning a TestValueGenerator of type T
     * \tparam T  type of values generated
     * \return TestValueGenerator<T>
     */
    template<typename T>
    static TestValueGenerator<T> testValueGenerator()
    {
        static const TestValueGenerator<T> testValueGenerator;
        return testValueGenerator;
    }

private:
    template<typename T>
    static const std::vector<T>& getTestVector();

    template<typename T>
    static std::enable_if_t<IsSerializableType<T>::value && !std::is_same<T, bool>::value, const T*>
    getBeginPointer();
    template<typename T>
    static std::enable_if_t<IsVectorOfSerializableType<T>::value, const T*> getBeginPointer();
    template<typename T>
    static std::enable_if_t<std::is_same<T, bool>::value, const T*> getBeginPointer();
    template<typename T>
    static std::enable_if_t<std::is_same<T, tensor>::value, const T*> getBeginPointer();

    template<typename T>
    static std::enable_if_t<IsSerializableType<T>::value && !std::is_same<T, bool>::value, const T*>
    getEndPointer();
    template<typename T>
    static std::enable_if_t<IsVectorOfSerializableType<T>::value, const T*> getEndPointer();
    template<typename T>
    static std::enable_if_t<std::is_same<T, bool>::value, const T*> getEndPointer();
    template<typename T>
    static std::enable_if_t<std::is_same<T, tensor>::value, const T*> getEndPointer();

    template<typename T>
    static std::enable_if_t<IsSerializableType<T>::value && !std::is_same<T, bool>::value, void>
    increment(const T** ptr);
    template<typename T>
    static std::enable_if_t<IsVectorOfSerializableType<T>::value, void> increment(const T** ptr);
    template<typename T>
    static std::enable_if_t<std::is_same<T, bool>::value, void> increment(const T** ptr);
    template<typename T>
    static std::enable_if_t<std::is_same<T, tensor>::value, void> increment(const T** ptr);

    static constexpr bool   testTrue    = true;
    static constexpr bool   testFalse   = false;
    static constexpr tensor testTensor1 = { { 1.6234, 2.4632, 3.1112 },
                                            { 4.66234, 5.9678, 6.088 },
                                            { 7.00001, 8.43535, 9.11233 } };
#if GMX_DOUBLE
    static constexpr tensor testTensor2 = { { 1, GMX_DOUBLE_EPS, 3 },
                                            { GMX_DOUBLE_MIN, 5, 6 },
                                            { 7, 8, GMX_DOUBLE_MAX } };
#else
    static constexpr tensor testTensor2 = { { 1, GMX_FLOAT_EPS, 3 },
                                            { GMX_FLOAT_MIN, 5, 6 },
                                            { 7, 8, GMX_FLOAT_MAX } };
#endif
};

// Begin implementations of TestValues methods
template<>
const std::vector<std::string>& TestValues::getTestVector()
{
    static const std::vector<std::string> testStrings{ "Test string\nwith newlines\n", "" };
    return testStrings;
}
template<>
const std::vector<int>& TestValues::getTestVector()
{
    static const std::vector<int> testInts{ 3, INT_MAX, INT_MIN };
    return testInts;
}
template<>
const std::vector<int64_t>& TestValues::getTestVector()
{
    static const std::vector<int64_t> testInt64s{ -7, LLONG_MAX, LLONG_MIN };
    return testInt64s;
}
template<>
const std::vector<float>& TestValues::getTestVector()
{
    static const std::vector<float> testFloats{ 33.9, GMX_FLOAT_MAX, GMX_FLOAT_MIN, GMX_FLOAT_EPS };
    return testFloats;
}
template<>
const std::vector<double>& TestValues::getTestVector()
{
    static const std::vector<double> testDoubles{ -123.45, GMX_DOUBLE_MAX, GMX_DOUBLE_MIN, GMX_DOUBLE_EPS };
    return testDoubles;
}

template<typename T>
std::enable_if_t<IsSerializableType<T>::value && !std::is_same<T, bool>::value, const T*> TestValues::getBeginPointer()
{
    return getTestVector<T>().data();
}
template<typename T>
std::enable_if_t<IsVectorOfSerializableType<T>::value, const T*> TestValues::getBeginPointer()
{
    return &getTestVector<typename T::value_type>();
}
template<typename T>
std::enable_if_t<std::is_same<T, bool>::value, const T*> TestValues::getBeginPointer()
{
    return &testTrue;
}
template<typename T>
std::enable_if_t<std::is_same<T, tensor>::value, const T*> TestValues::getBeginPointer()
{
    return &testTensor1;
}

template<typename T>
std::enable_if_t<IsSerializableType<T>::value && !std::is_same<T, bool>::value, const T*> TestValues::getEndPointer()
{
    return getTestVector<T>().data() + getTestVector<T>().size();
}
template<typename T>
std::enable_if_t<IsVectorOfSerializableType<T>::value, const T*> TestValues::getEndPointer()
{
    return &getTestVector<typename T::value_type>() + 1;
}
template<typename T>
std::enable_if_t<std::is_same<T, bool>::value, const T*> TestValues::getEndPointer()
{
    return nullptr;
}
template<typename T>
std::enable_if_t<std::is_same<T, tensor>::value, const T*> TestValues::getEndPointer()
{
    return nullptr;
}

template<typename T>
std::enable_if_t<IsSerializableType<T>::value && !std::is_same<T, bool>::value, void>
TestValues::increment(const T** ptr)
{
    ++(*ptr);
}
template<typename T>
std::enable_if_t<IsVectorOfSerializableType<T>::value, void> TestValues::increment(const T** ptr)
{
    ++(*ptr);
}
template<typename T>
std::enable_if_t<std::is_same<T, bool>::value, void> TestValues::increment(const T** ptr)
{
    *ptr = (*ptr == &testTrue) ? &testFalse : nullptr;
}
template<typename T>
std::enable_if_t<std::is_same<T, tensor>::value, void> TestValues::increment(const T** ptr)
{
    *ptr = (*ptr == &testTensor1) ? &testTensor2 : nullptr;
}

template<typename T>
typename TestValues::TestValueGenerator<T>::Iterator TestValues::TestValueGenerator<T>::begin() const
{
    return TestValues::TestValueGenerator<T>::Iterator(getBeginPointer<T>());
}

template<typename T>
typename TestValues::TestValueGenerator<T>::Iterator TestValues::TestValueGenerator<T>::end() const
{
    return TestValues::TestValueGenerator<T>::Iterator(getEndPointer<T>());
}

template<typename T>
typename TestValues::TestValueGenerator<T>::Iterator TestValues::TestValueGenerator<T>::Iterator::operator++()
{
    TestValues::increment(&ptr_);
    return *this;
}
// End implementations of TestValues methods

//! Write scalar input to CheckpointData
template<typename T>
typename std::enable_if_t<IsSerializableType<T>::value, void>
writeInput(const std::string& key, const T& inputValue, WriteCheckpointData* checkpointData)
{
    checkpointData->scalar(key, &inputValue);
}
//! Read scalar from CheckpointData and test if equal to input
template<typename T>
typename std::enable_if_t<IsSerializableType<T>::value, void>
testOutput(const std::string& key, const T& inputValue, ReadCheckpointData* checkpointData)
{
    T outputValue;
    checkpointData->scalar(key, &outputValue);
    EXPECT_EQ(inputValue, outputValue);
}
//! Write vector input to CheckpointData
template<typename T>
void writeInput(const std::string& key, const std::vector<T>& inputVector, WriteCheckpointData* checkpointData)
{
    checkpointData->arrayRef(key, makeConstArrayRef(inputVector));
}
//! Read vector from CheckpointData and test if equal to input
template<typename T>
void testOutput(const std::string& key, const std::vector<T>& inputVector, ReadCheckpointData* checkpointData)
{
    std::vector<T> outputVector;
    outputVector.resize(inputVector.size());
    checkpointData->arrayRef(key, makeArrayRef(outputVector));
    EXPECT_THAT(outputVector, ::testing::ContainerEq(inputVector));
}
//! Write tensor input to CheckpointData
void writeInput(const std::string& key, const tensor inputTensor, WriteCheckpointData* checkpointData)
{
    checkpointData->tensor(key, inputTensor);
}
//! Read tensor from CheckpointData and test if equal to input
void testOutput(const std::string& key, const tensor inputTensor, ReadCheckpointData* checkpointData)
{
    tensor outputTensor = { { 0 } };
    checkpointData->tensor(key, outputTensor);
    std::array<std::array<real, 3>, 3> inputTensorArray = {
        { { inputTensor[XX][XX], inputTensor[XX][YY], inputTensor[XX][ZZ] },
          { inputTensor[YY][XX], inputTensor[YY][YY], inputTensor[YY][ZZ] },
          { inputTensor[ZZ][XX], inputTensor[ZZ][YY], inputTensor[ZZ][ZZ] } }
    };
    std::array<std::array<real, 3>, 3> outputTensorArray = {
        { { outputTensor[XX][XX], outputTensor[XX][YY], outputTensor[XX][ZZ] },
          { outputTensor[YY][XX], outputTensor[YY][YY], outputTensor[YY][ZZ] },
          { outputTensor[ZZ][XX], outputTensor[ZZ][YY], outputTensor[ZZ][ZZ] } }
    };
    EXPECT_THAT(outputTensorArray, ::testing::ContainerEq(inputTensorArray));
}

/*! \internal
 * \brief CheckpointData test fixture
 *
 * Test whether input is equal to output, either with a single data type
 * or with a combination of three data types.
 */
class CheckpointDataTest : public ::testing::Test
{
public:
    using WriteFunction = std::function<void(WriteCheckpointData*)>;
    using TestFunction  = std::function<void(ReadCheckpointData*)>;

    // List of functions to write values to checkpoint
    std::vector<WriteFunction> writeFunctions_;
    // List of functions to test read checkpoint object
    std::vector<TestFunction> testFunctions_;

    // Add values to write / test lists
    template<typename T>
    void addTestValues()
    {
        for (const auto& inputValue : TestValues::testValueGenerator<T>())
        {
            std::string key = "value" + std::to_string(writeFunctions_.size());
            writeFunctions_.emplace_back([key, inputValue](WriteCheckpointData* checkpointData) {
                writeInput(key, inputValue, checkpointData);
            });
            testFunctions_.emplace_back([key, inputValue](ReadCheckpointData* checkpointData) {
                testOutput(key, inputValue, checkpointData);
            });
        }
    }

    /* This shuffles the write and test functions (as writing and reading can happen
     * if different orders), then writes all data to a CheckpointData object.
     * The CheckpointData object is serialized to file, and then re-read into
     * a new CheckpointData object. The test functions are then used to assert
     * that all data is present in the new object.
     */
    void test()
    {
        /* Randomize order of writing and testing - checkpoint data can be
         * accessed in any order!
         * Having the same order makes this reproducible, so at least for now we're
         * ok seeding the rng with default value and silencing clang-tidy */
        // NOLINTNEXTLINE(cert-msc51-cpp)
        auto rng = std::default_random_engine{};
        std::shuffle(std::begin(writeFunctions_), std::end(writeFunctions_), rng);
        std::shuffle(std::begin(testFunctions_), std::end(testFunctions_), rng);

        // Write value to CheckpointData & serialize
        {
            WriteCheckpointDataHolder writeCheckpointDataHolder;
            auto writeCheckpointData = writeCheckpointDataHolder.checkpointData("test");
            for (const auto& writeFunction : writeFunctions_)
            {
                writeFunction(&writeCheckpointData);
            }

            auto*               file = gmx_fio_open(filename_, "w");
            FileIOXdrSerializer serializer(file);
            writeCheckpointDataHolder.serialize(&serializer);
            gmx_fio_close(file);
        }

        // Deserialize values and test against reference
        {
            auto*               file = gmx_fio_open(filename_, "r");
            FileIOXdrSerializer deserializer(file);

            ReadCheckpointDataHolder readCheckpointDataHolder;
            readCheckpointDataHolder.deserialize(&deserializer);
            gmx_fio_close(file);

            auto readCheckpointData = readCheckpointDataHolder.checkpointData("test");
            for (const auto& testFunction : testFunctions_)
            {
                testFunction(&readCheckpointData);
            }
        }

        // Object can be reused
        writeFunctions_.clear();
        testFunctions_.clear();
    }

    // The different functions to add test values - in a list to simplify looping over them
    std::vector<std::function<void()>> addTestValueFunctions_ = {
        [this]() { addTestValues<std::string>(); },
        [this]() { addTestValues<int>(); },
        [this]() { addTestValues<int64_t>(); },
        [this]() { addTestValues<bool>(); },
        [this]() { addTestValues<float>(); },
        [this]() { addTestValues<double>(); },
        [this]() { addTestValues<std::vector<std::string>>(); },
        [this]() { addTestValues<std::vector<int>>(); },
        [this]() { addTestValues<std::vector<int64_t>>(); },
        [this]() { addTestValues<std::vector<float>>(); },
        [this]() { addTestValues<std::vector<double>>(); },
        [this]() { addTestValues<tensor>(); }
    };

    // The types we're testing - for scoped trace output only
    std::vector<std::string> testingTypes_ = { "std::string",
                                               "int",
                                               "int64_t",
                                               "bool",
                                               "float",
                                               "double",
                                               "std::vector<std::string>",
                                               "std::vector<int>",
                                               "std::vector<int64_t>",
                                               "std::vector<float>",
                                               "std::vector<double>",
                                               "tensor" };

    // We'll need a temporary file to write / read our dummy checkpoint to
    TestFileManager       fileManager_;
    std::filesystem::path filename_ = fileManager_.getTemporaryFilePath("test.cpt");
};

TEST_F(CheckpointDataTest, SingleDataTest)
{
    // Test data types separately
    const int numTypes = addTestValueFunctions_.size();
    for (int type = 0; type < numTypes; ++type)
    {
        SCOPED_TRACE("Using test values of type " + testingTypes_[type]);
        addTestValueFunctions_[type]();
        test();
    }
}

TEST_F(CheckpointDataTest, MultiDataTest)
{
    // All combinations of 2 different data types
    const int numTypes = addTestValueFunctions_.size();
    for (int type1 = 0; type1 < numTypes; ++type1)
    {
        for (int type2 = type1; type2 < numTypes; ++type2)
        {
            SCOPED_TRACE("Using test values of type " + testingTypes_[type1] + " and "
                         + testingTypes_[type2]);
            addTestValueFunctions_[type1]();
            addTestValueFunctions_[type2]();
            test();
        }
    }
}

} // namespace
} // namespace gmx::test
