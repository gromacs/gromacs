/*
 * Copyright © 2010-2013  Felix Höfling
 * Copyright © 2013-2014  Manuel Dibak
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#define BOOST_TEST_MODULE h5xx_attribute
#include <boost/test/unit_test.hpp>
#include <boost/array.hpp>

#include <h5xx/attribute.hpp>
#include <h5xx/group.hpp>
#include <h5xx/policy/string.hpp>

#include <test/ctest_full_output.hpp>
#include <test/catch_boost_no_throw.hpp>
#include <test/fixture.hpp>

#include <stdio.h>

BOOST_GLOBAL_FIXTURE( ctest_full_output );

namespace fixture { // preferred over BOOST_FIXTURE_TEST_SUITE

char filename[] = "test_h5xx_attribute.h5";
typedef h5file<filename> BOOST_AUTO_TEST_CASE_FIXTURE;

using namespace h5xx;

BOOST_AUTO_TEST_CASE( construction )
{
    BOOST_CHECK_NO_THROW(attribute());                 // default constructor

    write_attribute(file, "foo", 1);                   // create attribute in a file's root group
    BOOST_CHECK_NO_THROW(attribute(file, "foo"));      // open existing attribute on-the-fly

    attribute foo(file, "foo");
//    BOOST_CHECK_EQUAL(get_name(foo), "/");           // path of the object the attribute is attached to
//    BOOST_CHECK_EQUAL(foo.name(), "foo");            // name of the attribute
    BOOST_CHECK_EQUAL(get_name(foo), "/foo");          // NEW: full path of the attribute
    BOOST_CHECK_EQUAL(foo.name(), "/foo");             // NEW: full path of the attribute
    BOOST_CHECK(foo.valid());

    hid_t hid = foo.hid();
    attribute bar;
    BOOST_CHECK_THROW(bar = foo, h5xx::error);         // assignment is not allowed (it makes a copy)
    BOOST_CHECK_NO_THROW(bar = move(foo));             // move assignment
    BOOST_CHECK_EQUAL(bar.hid(), hid);
    BOOST_CHECK(!foo.valid());

    BOOST_CHECK_THROW(attribute g(bar), h5xx::error);  // copying is not allowed
    BOOST_CHECK_NO_THROW(attribute g(move(bar)));      // copying from temporary is allowed (with move semantics)
    BOOST_CHECK(!bar.valid());
}

BOOST_AUTO_TEST_CASE( scalar_fundamental )
{

    bool test_bool = true;
    BOOST_CHECK_NO_THROW(write_attribute(file, "bool, scalar", test_bool));
    BOOST_CHECK(read_attribute<bool>(file, "bool, scalar") == test_bool);
    BOOST_CHECK(exists_attribute(file, "bool, scalar") == true);
    BOOST_CHECK_NO_THROW(delete_attribute(file, "bool, scalar"));
    BOOST_CHECK(exists_attribute(file, "bool, scalar") == false);
    H5E_BEGIN_TRY{
        BOOST_CHECK_THROW(read_attribute<bool>(file, "bool, scalar"), h5xx::error);
    } H5E_END_TRY

    double double_value = std::sqrt(2.L);
    BOOST_CHECK_NO_THROW(write_attribute(file, "double, scalar", 2));   // store something of wrong type first
    BOOST_CHECK_NO_THROW(write_attribute(file, "double, scalar", double_value));
    BOOST_CHECK(exists_attribute(file, "double, scalar") == true);
    BOOST_CHECK(read_attribute<double>(file,"double, scalar") == double_value);

    uint64_t uint_value = 9223372036854775783LLU;  // largest prime below 2^63
    BOOST_CHECK_NO_THROW(write_attribute(file, "integral, scalar", uint64_t(1)));   // store wrong value of correct type first
    BOOST_CHECK_NO_THROW(write_attribute(file, "integral, scalar", uint_value));  // overwrite value
    BOOST_CHECK(read_attribute<uint64_t>(file, "integral, scalar") == uint_value);

}

BOOST_AUTO_TEST_CASE( scalar_stdstring )
{
    std::string stdstring = "HALMD";
    std::string wrongstring = "wrong";
    std::string read;
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, stdstring", wrongstring));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, stdstring", stdstring));
    BOOST_CHECK(exists_attribute(file, "scalar, stdstring"));
    BOOST_CHECK_NO_THROW(read_attribute<std::string>(file, "scalar, stdstring"));
    BOOST_CHECK_NO_THROW(read = read_attribute<std::string>(file, "scalar, stdstring"));
    BOOST_CHECK_EQUAL(read, stdstring);
    BOOST_CHECK_EQUAL(read.compare(0, 5, stdstring), 0);
    BOOST_CHECK_EQUAL(read.size(), stdstring.size());
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, stdstring, nullterm", stdstring, policy::string::null_terminated()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, stdstring, nullpad", stdstring, policy::string::null_padded()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, stdstring, spacepad", stdstring, policy::string::space_padded()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, stdstring, varlen", stdstring, policy::string::variable_length()));
    BOOST_CHECK(stdstring == read_attribute<std::string>(file, "scalar, stdstring, varlen"));
    BOOST_CHECK(stdstring == read_attribute<std::string>(file, "scalar, stdstring, nullpad"));
    BOOST_CHECK(stdstring == read_attribute<std::string>(file, "scalar, stdstring, spacepad"));
    BOOST_CHECK(stdstring == read_attribute<std::string>(file, "scalar, stdstring, nullterm"));
/*    typedef hid_t (*write_function_type)(h5xx::file const&, std::string const&, std::string const&);
    typedef hid_t (*write_function_type2)(h5xx::file const&, std::string const&, std::string const&, policy::string::null_terminated);
    write_function_type f1 = &write_attribute<std::string, h5xx::file>;
    write_function_type f2 = &write_attribute<std::string, h5xx::file, policy::string::null_terminated>;
//     BOOST_CHECK_EQUAL(f1, f2);
    */
}

BOOST_AUTO_TEST_CASE( scalar_cstring )
{
    char const* cstring = "Highly accelerated large-scale molecular dynamics simulation package";;
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring", "wrong string"));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring", cstring));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring, nullterm", cstring));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring, varlen", cstring, policy::string::variable_length()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring, nullpad", cstring, policy::string::null_padded()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring, spacepad", cstring, policy::string::space_padded()));
//    BOOST_CHECK_NO_THROW(write_attribute(file, "scalar, cstring, empty", ""));

    BOOST_CHECK(exists_attribute(file, "scalar, cstring"));
    BOOST_CHECK(read_attribute<std::string>(file, "scalar, cstring") == cstring);
    BOOST_CHECK(read_attribute<std::string>(file, "scalar, cstring, nullterm") == cstring);
    BOOST_CHECK(read_attribute<std::string>(file, "scalar, cstring, varlen") == cstring);
    BOOST_CHECK(read_attribute<std::string>(file, "scalar, cstring, nullpad") == cstring);
    BOOST_CHECK(read_attribute<std::string>(file, "scalar, cstring, spacepad") == cstring);
//    BOOST_CHECK(read_attribute<std::string>(file, "scalar, cstring, empty") == "");
}

BOOST_AUTO_TEST_CASE( array_int )
{
    typedef boost::array<int, 2> int_array_type;
    typedef boost::array<double, 2> double_array_type;
    double_array_type wrong_array = {{1, 2}};
    int_array_type int_array = {{ 1233, 12344}};
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, int", wrong_array));          // write wrong array first
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, int", int_array));            // overwrite with correct array
    BOOST_CHECK(exists_attribute(file, "array, int"));
    BOOST_CHECK(read_attribute<int_array_type>(file, "array, int") == int_array);    // check if overwriting was successful
}

BOOST_AUTO_TEST_CASE( array_bool )
{
    typedef boost::array<bool, 2> bool_array_type;
    bool_array_type true_array = {{true, false}};
    bool_array_type false_array = {{false, true}};
    bool_array_type read;
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, bool", false_array ));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, bool", true_array));
    BOOST_CHECK_NO_THROW(read = read_attribute<bool_array_type>(file, "array, bool"));
    BOOST_CHECK(read == true_array);
}

BOOST_AUTO_TEST_CASE( array_cstring )
{
    typedef boost::array<char const*, 3> cstring_array_type;
    typedef boost::array<std::string, 3> string_array_type;
    cstring_array_type cstring_array = {{
        "HAL's MD package",
        "Highly accelerated large-scale molecular dynamics simulation package",
        "HALMD"
    }};
    string_array_type string_array = {{
        "HAL's MD package",
        "Highly accelerated large-scale molecular dynamics simulation package",
        "HALMD"
    }};
    string_array_type read;
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, cstring", cstring_array));
    BOOST_CHECK_NO_THROW(read = read_attribute<string_array_type>(file, "array, cstring"));
    BOOST_CHECK(read == string_array);
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, cstring, nullterm", cstring_array, policy::string::null_terminated()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, cstring, varlen", cstring_array, policy::string::variable_length()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, cstring, nullpad", cstring_array, policy::string::null_padded()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, cstring, spacepad", cstring_array, policy::string::space_padded()));
    BOOST_CHECK(read_attribute<string_array_type>(file, "array, cstring, nullterm") == string_array);
    BOOST_CHECK(read_attribute<string_array_type>(file, "array, cstring, varlen") == string_array);
    BOOST_CHECK(read_attribute<string_array_type>(file, "array, cstring, nullpad") == string_array);
    BOOST_CHECK(read_attribute<string_array_type>(file, "array, cstring, spacepad") == string_array);
}

BOOST_AUTO_TEST_CASE( array_string )
{
    typedef boost::array<std::string, 3> string_array_type;
    string_array_type string_array = {{
        "HAL's MD package",
        "Highly accelerated large-scale molecular dynamics simulation package",
        "HALMD"
    }};
    string_array_type read;
    string_array_type read_nullterm;
    string_array_type read_variable_length;
    string_array_type read_nullpad;
    string_array_type read_spacepad;
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, string", string_array));
    BOOST_CHECK_NO_THROW(read = read_attribute<string_array_type>(file, "array, string"));
    BOOST_CHECK(read == string_array);
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, string, nullterm", string_array, policy::string::null_terminated()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, string, nullpad", string_array, policy::string::null_padded()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, string, spacepad", string_array, policy::string::space_padded()));
    BOOST_CHECK_NO_THROW(write_attribute(file, "array, string, varlen", string_array, policy::string::variable_length()));
    BOOST_CHECK_NO_THROW(read_nullterm = read_attribute<string_array_type>(file, "array, string, nullterm"));
    BOOST_CHECK_NO_THROW(read_variable_length = read_attribute<string_array_type>(file, "array, string, varlen"));
    BOOST_CHECK_NO_THROW(read_nullpad = read_attribute<string_array_type>(file, "array, string, nullpad"));
    BOOST_CHECK_NO_THROW(read_spacepad = read_attribute<string_array_type>(file, "array, string, spacepad"));
    BOOST_CHECK(string_array == read_nullterm);
    BOOST_CHECK(string_array == read_nullpad);
    BOOST_CHECK(string_array == read_spacepad);
    BOOST_CHECK(string_array == read_variable_length);
}

BOOST_AUTO_TEST_CASE( stdvector_double )
{
    typedef std::vector<double> double_vector_type;
    static const double values[] = {std::sqrt(2), 3, 77};
    double_vector_type double_vector(values,values + sizeof(values)/sizeof(double));
    double_vector_type read;
    BOOST_CHECK_NO_THROW(write_attribute(file, "stdvector, double", double_vector));
    BOOST_CHECK_NO_THROW(read = read_attribute<double_vector_type>(file, "stdvector, double"));
    BOOST_CHECK(read == double_vector);
}

BOOST_AUTO_TEST_CASE( stdvector_int )
{
    typedef std::vector<int> int_vector_type;
    static const int values[] = {1, 2, 3};
    int_vector_type int_vector(values, values + sizeof(values)/sizeof(int));
    int_vector_type read;
    BOOST_CHECK_NO_THROW(write_attribute(file, "stdvector, int", int_vector));
    BOOST_CHECK_NO_THROW(read = read_attribute<int_vector_type>(file, "stdvector, int"));
    BOOST_CHECK(read == int_vector);
}

BOOST_AUTO_TEST_CASE( stdvector_string )
{
    std::vector<std::string> output, input;
    output.push_back("HAL's MD package");
    output.push_back("Highly accelerated large-scale molecular dynamics simulation package");
    output.push_back("HALMD");
    BOOST_CHECK_NO_THROW(write_attribute(file, "std::vector, string", output));
    BOOST_CHECK_NO_THROW(input = read_attribute<std::vector<std::string> >(file, "std::vector, string"));
    BOOST_CHECK_EQUAL_COLLECTIONS(input.begin(), input.end(), output.begin(), output.end());
}

BOOST_AUTO_TEST_CASE( boost_multi_array)
{
    typedef boost::multi_array<int, 3> multi_array3;
    int data3[] = {99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76};
    multi_array3 multi_array_value(boost::extents[2][3][4]);
    multi_array_value.assign(data3, data3 + 2 * 3 * 4);
    multi_array3 read(boost::extents[2][3][4]);
    BOOST_CHECK_NO_THROW(write_attribute(file, "boost multi array, int", multi_array_value));
    BOOST_CHECK_NO_THROW(read = read_attribute<multi_array3>(file, "boost multi array, int"));
    BOOST_CHECK(read == multi_array_value);
}
} //namespace fixture
