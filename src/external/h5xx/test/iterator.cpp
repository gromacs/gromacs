/*
 * Copyright © 2018      Matthias Werner
 * Copyright © 2018-2019 Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

/**
 * Requirements of Iterator(It)
 * http://en.cppreference.com/w/cpp/concept/ForwardIterator
 *
 *  - DefaultConstructible
 *  - CopyConstructible
 *  - CopyAssignable
 *  - EqualityComparable
 *  - Destructible
 *  - lvalues are Swappable
 *  - std::iterator_traits<It> has member typedefs:
 *    value_type, difference_type, reference, pointer and iterator_category
 *
 *  and
 *
 *  - operator* is defined - tested
 *  - operator++() is defined (returns reference to iterator) - tested
 *  - operator++(int) - tested
 *  - (void) i++ equivalent (void)++i - tested
 *  - operator!= / operator== - tested
 *  - It->m should be equivalent to (*It).m -tested
 *
 *  and (for ForwardIterator)
 *
 *  - Multipass guarantee - tested
 */
#define BOOST_TEST_MODULE h5xx_group
#include <boost/test/unit_test.hpp>

#include <h5xx/group.hpp>
#include <h5xx/dataset.hpp>

#include <test/ctest_full_output.hpp>
#include <test/catch_boost_no_throw.hpp>
#include <test/fixture.hpp>

BOOST_GLOBAL_FIXTURE( ctest_full_output );

namespace fixture { // preferred over BOOST_FIXTURE_TEST_SUITE

char filename[] = "test_h5xx_iterator.h5";
typedef h5file<filename> BOOST_AUTO_TEST_CASE_FIXTURE;

using namespace h5xx;

BOOST_AUTO_TEST_CASE( iterator_constructors )
{
    BOOST_TEST_MESSAGE("\nTesting iterator constructors:");

    group container_group(file);

    // default constructor
    BOOST_TEST_MESSAGE("  default constructor");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator dset_iter);
    BOOST_CHECK_NO_THROW(container<group>::iterator sgroup_iter);

    container<dataset>::iterator dset_iter;
    container<group>::iterator sgroup_iter;

    // custom constructor
    BOOST_TEST_MESSAGE("  constructor on empty group");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator container_group);
    BOOST_CHECK_NO_THROW(container<group>::iterator container_group);

    // populate group
    create_dataset<int>(container_group, "dset");
    group(container_group, "grp");

    // custom constructor on non-empty group
    BOOST_TEST_MESSAGE("  constructor on non-empty group");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator container_group);
    BOOST_CHECK_NO_THROW(container<group>::iterator container_group);

    // copy constructor
    BOOST_TEST_MESSAGE("  copy constructor");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator dset_iter_2 = dset_iter);
    BOOST_CHECK_NO_THROW(container<group>::iterator sgroup_iter_2 = sgroup_iter);

    // move constructor
    BOOST_TEST_MESSAGE("  move constructor");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator dset_iter_2 = container_group.datasets().begin());    // temporary on the right
    BOOST_CHECK_NO_THROW(container<group>::iterator sgroup_iter_2 = container_group.groups().begin());

    container<dataset>::iterator dset_iter_2(container_group);
    container<group>::iterator sgroup_iter_2(container_group);

    // assignment operator
    BOOST_TEST_MESSAGE("  assignment operator");
    BOOST_CHECK_NO_THROW(dset_iter = dset_iter_2);
    BOOST_CHECK_NO_THROW(sgroup_iter_2 = sgroup_iter);
    BOOST_CHECK_NO_THROW(dset_iter = container<dataset>::iterator());
    BOOST_CHECK_NO_THROW(sgroup_iter = container<group>::iterator(container_group));
}

BOOST_AUTO_TEST_CASE( iterator_requirements )
{
    BOOST_TEST_MESSAGE("\nTesting basic requirements of ForwardIterator:");

    // DefaultConstructible
    BOOST_TEST_MESSAGE("  DefaultConstructible");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator());
    BOOST_CHECK_NO_THROW(container<group>::iterator());

    // MoveConstructible
    BOOST_TEST_MESSAGE("  MoveConstructible");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator temp{container<dataset>::iterator()});
    BOOST_CHECK_NO_THROW(container<group>::iterator temp{container<group>::iterator()});

    container<dataset>::iterator dset_iter, dset_iter_2, dset_iter_3;
    container<group>::iterator sgroup_iter, sgroup_iter_2, sgroup_iter_3;

    // MoveAssignable
    BOOST_TEST_MESSAGE("  MoveAssignable");
    BOOST_CHECK_NO_THROW(dset_iter = container<dataset>::iterator());
    BOOST_CHECK_NO_THROW(sgroup_iter = container<group>::iterator());

    // CopyConstructible
    BOOST_TEST_MESSAGE("  CopyConstructible");
    BOOST_CHECK_NO_THROW(container<dataset>::iterator dset_iter);
    BOOST_CHECK_NO_THROW(container<group>::iterator sgroup_iter);

    // CopyAssignable
    BOOST_TEST_MESSAGE("  CopyAssignable");
    BOOST_CHECK_NO_THROW(dset_iter_2 = dset_iter;);
    BOOST_CHECK_NO_THROW(sgroup_iter_2 = sgroup_iter);

    // EqualityComparable
    BOOST_TEST_MESSAGE("  EqualityComparable");
    BOOST_CHECK(dset_iter == dset_iter_2);
    BOOST_CHECK(sgroup_iter == sgroup_iter_2);

    // lvalues are Swappable
    BOOST_TEST_MESSAGE("  lvalues are swappable");
    BOOST_CHECK_NO_THROW(std::swap(dset_iter, dset_iter_2));
    BOOST_CHECK_NO_THROW(std::swap(sgroup_iter, sgroup_iter_2));

    // member typedefs value_type, difference_type, reference, pointer and iterator_category
    BOOST_TEST_MESSAGE("  std::iterator_traits for non-const iterator");
    BOOST_CHECK(typeid(container<dataset>::iterator::value_type) == typeid(dataset));
    BOOST_CHECK(typeid(container<group>::iterator::value_type) == typeid(group));

    BOOST_CHECK(typeid(container<dataset>::iterator::difference_type) == typeid(std::ptrdiff_t));
    BOOST_CHECK(typeid(container<group>::iterator::difference_type) ==  typeid(std::ptrdiff_t));

    BOOST_CHECK(typeid(container<dataset>::iterator::reference) == typeid(dataset&));
    BOOST_CHECK(typeid(container<group>::iterator::reference) == typeid(group&));

    BOOST_CHECK(typeid(container<dataset>::iterator::pointer) == typeid(dataset*));
    BOOST_CHECK(typeid(container<group>::iterator::pointer) == typeid(group*));

    BOOST_CHECK(typeid(container<dataset>::iterator::iterator_category) == typeid(std::forward_iterator_tag));
    BOOST_CHECK(typeid(container<group>::iterator::iterator_category) == typeid(std::forward_iterator_tag));

    BOOST_TEST_MESSAGE("  std::iterator_traits for const iterator");
    BOOST_CHECK(typeid(container<dataset>::const_iterator::value_type) == typeid(dataset));
    BOOST_CHECK(typeid(container<group>::const_iterator::value_type) == typeid(group));

    BOOST_CHECK(typeid(container<dataset>::const_iterator::difference_type) == typeid(std::ptrdiff_t));
    BOOST_CHECK(typeid(container<group>::const_iterator::difference_type) ==  typeid(std::ptrdiff_t));

    BOOST_CHECK(typeid(container<dataset>::const_iterator::reference) == typeid(dataset const&));
    BOOST_CHECK(typeid(container<group>::const_iterator::reference) == typeid(group const&));

    BOOST_CHECK(typeid(container<dataset>::const_iterator::pointer) == typeid(dataset const*));
    BOOST_CHECK(typeid(container<group>::const_iterator::pointer) == typeid(group const*));

    BOOST_CHECK(typeid(container<dataset>::const_iterator::iterator_category) == typeid(std::forward_iterator_tag));
    BOOST_CHECK(typeid(container<group>::const_iterator::iterator_category) == typeid(std::forward_iterator_tag));

    // further test with a real group
    group container_group(file);
    create_dataset<int>(container_group, "dset1");
    create_dataset<int>(container_group, "dset2");
    group(container_group, "grp1");
    group(container_group, "grp2");
    group(container_group, "grp3");

    container<group>::iterator sgroup_multipass_1 = container_group.groups().begin();
    container<group>::iterator sgroup_multipass_2 = container_group.groups().begin();

    // Multipass guarantee
    BOOST_TEST_MESSAGE("  Multipass guarantee");

    auto datasets = container_group.datasets();
    container<dataset>::iterator dset_it1 = datasets.begin();
    for (auto dset_it2 = datasets.begin(); dset_it2 != datasets.end(); ++dset_it2) {
        BOOST_CHECK(dset_it1 == dset_it2);
        BOOST_CHECK(dset_it1.get_name() == dset_it2.get_name());
        dset_it1++;
    }

    auto groups = container_group.groups();
    container<group>::iterator grp_it1 = groups.begin();
    for (auto grp_it2 = groups.begin(); grp_it2 != groups.end(); ++grp_it2) {
        BOOST_CHECK(grp_it1 == grp_it2);
        BOOST_CHECK(grp_it1.get_name() == grp_it2.get_name());
        grp_it1++;
    }
}

BOOST_AUTO_TEST_CASE( iterator_expressions )
{
    BOOST_TEST_MESSAGE("\nTesting basic iterator expressions:");

    // create a small group with a single dataset and no subgroups
    group container_group(file);
    create_dataset<int>(container_group, "dset");
    create_dataset<double>(container_group, "dset2");

    // a default constructed iterator, and one over the above group
    container<dataset>::iterator dset_iter, dset_iter_2(container_group);
    container<group>::iterator sgroup_iter, sgroup_iter_2(container_group);

    // operator*
    BOOST_TEST_MESSAGE("  dereference operator");
    // throws h5xx::error since default constructed (underlying group does not exist)
    BOOST_CHECK_THROW(*sgroup_iter, h5xx::error);
    BOOST_CHECK_THROW(*dset_iter, h5xx::error);
    BOOST_CHECK_EQUAL(get_name(*dset_iter_2), "/dset");
    BOOST_CHECK_THROW(*sgroup_iter_2, std::out_of_range);

    // operator++()
    BOOST_TEST_MESSAGE("  pre-increment operator");
    BOOST_CHECK_THROW(++dset_iter, h5xx::error);
    BOOST_CHECK_THROW(++sgroup_iter, h5xx::error);
    BOOST_CHECK_EQUAL(get_name(*(++dset_iter_2)), "/dset2");
    BOOST_CHECK_THROW(*(++sgroup_iter_2), std::out_of_range);

    auto datasets = container_group.datasets();
    auto groups = container_group.groups();
    dset_iter_2 = datasets.begin();                 // go back to start
    sgroup_iter_2 = groups.begin();

    // operator++(int)
    BOOST_TEST_MESSAGE("  post-increment operator");
    BOOST_CHECK_THROW(dset_iter++, h5xx::error);
    BOOST_CHECK_THROW(sgroup_iter++, h5xx::error);
    BOOST_CHECK_EQUAL(get_name(*dset_iter_2++), "/dset");
    BOOST_CHECK_THROW(*sgroup_iter_2++, std::out_of_range);

    // operator!= and operator==
    BOOST_TEST_MESSAGE("  comparison operators");
    BOOST_CHECK_THROW(dset_iter == dset_iter_2, h5xx::error);   // different parent groups
    BOOST_CHECK(groups.begin() == groups.end());                // empty group
    BOOST_CHECK(groups.begin() == groups.cend());               // non-const and const

    dset_iter_2 = datasets.begin();                             // go back to start
    BOOST_CHECK(dset_iter_2 == datasets.begin());               // dset_iter_2 was never accessed
    BOOST_CHECK(dset_iter_2 == datasets.begin());               // dset_iter_2 accessed before for comparison

    dset_iter_2 = datasets.begin();                             // reset iterator and dereference
    BOOST_CHECK_NO_THROW(*dset_iter_2);
    BOOST_CHECK(datasets.begin() == dset_iter_2);               // dset_iter_2 dereferenced before
    BOOST_CHECK(++dset_iter_2 != datasets.end());

    // (void) i++ == (void) ++i
    auto dset_iter_3 = datasets.begin(); dset_iter_3++;         // post-increment here, pre-increment there
    BOOST_CHECK(dset_iter_2 == dset_iter_3);
    BOOST_CHECK_EQUAL(get_name(*dset_iter_2), get_name(*dset_iter_3));

    // check It->m equivalent to (*It).m
    BOOST_TEST_MESSAGE("  equivalence of operator-> and operator*");
    BOOST_CHECK_EQUAL((*dset_iter_2).hid(), dset_iter_2->hid());
}

BOOST_AUTO_TEST_CASE( default_group )
{
    BOOST_TEST_MESSAGE("\nTesting iterator over default constructed group");
    group container_group;
    container<dataset> datasets = container_group.datasets();
    container<group> groups = container_group.groups();

    BOOST_TEST_MESSAGE("  setting iterators to begin of containers");
    container<dataset>::iterator dset_iter;
    container<dataset>::const_iterator dset_citer;
    container<group>::iterator sgroup_iter;
    container<group>::const_iterator sgroup_citer;

    BOOST_CHECK_NO_THROW(dset_iter = datasets.begin());
    BOOST_CHECK_NO_THROW(dset_citer = datasets.cbegin());
    BOOST_CHECK_NO_THROW(sgroup_iter = groups.begin());
    BOOST_CHECK_NO_THROW(sgroup_citer = groups.cbegin());

    BOOST_TEST_MESSAGE("  equality of begin/cbegin and end/cend iterators");
    BOOST_CHECK(dset_iter == datasets.end());
    BOOST_CHECK(dset_citer == datasets.cend());
    BOOST_CHECK(sgroup_iter == groups.end());
    BOOST_CHECK(sgroup_citer == groups.cend());

    BOOST_CHECK_THROW(*dset_iter, std::out_of_range);
    BOOST_CHECK_THROW(*dset_citer, std::out_of_range);
    BOOST_CHECK_THROW(*sgroup_iter, std::out_of_range);
    BOOST_CHECK_THROW(*sgroup_citer, std::out_of_range);
}

BOOST_AUTO_TEST_CASE( empty_group )
{
    BOOST_TEST_MESSAGE("\nTesting iterator over empty group");
    group container_group(file);
    container<dataset> datasets = container_group.datasets();
    container<group> groups = container_group.groups();

    BOOST_TEST_MESSAGE("  setting iterators to begin of containers");
    container<dataset>::iterator dset_iter;
    container<dataset>::const_iterator dset_citer;
    container<group>::iterator sgroup_iter;
    container<group>::const_iterator sgroup_citer;

    BOOST_CHECK_NO_THROW(dset_iter = datasets.begin());
    BOOST_CHECK_NO_THROW(dset_citer = datasets.cbegin());
    BOOST_CHECK_NO_THROW(sgroup_iter = groups.begin());
    BOOST_CHECK_NO_THROW(sgroup_citer = groups.cbegin());

    // begin and end iterators should be equal in empty group
    BOOST_TEST_MESSAGE("  equality of begin/cbegin and end/cend iterators");
    BOOST_CHECK(dset_iter == datasets.end());
    BOOST_CHECK(dset_citer == datasets.cend());
    BOOST_CHECK(sgroup_iter == groups.end());
    BOOST_CHECK(sgroup_citer == groups.cend());

    // dereferencing of past-the-end iterators is not allowed
    BOOST_CHECK_THROW(*dset_iter, std::out_of_range);
    BOOST_CHECK_THROW(*sgroup_iter, std::out_of_range);
    BOOST_CHECK_THROW(*dset_citer, std::out_of_range);
    BOOST_CHECK_THROW(*sgroup_citer, std::out_of_range);
}

BOOST_AUTO_TEST_CASE( only_datasets )
{
    BOOST_TEST_MESSAGE("\nTesting iterator over dataset-only group");
    group container_group(file);

    create_dataset<int>(container_group, "dset1");
    create_dataset<int>(container_group, "dset2");

    container<dataset> datasets = container_group.datasets();
    container<group> groups = container_group.groups();

    // add one more after constructing the container adapter
    create_dataset<double>(container_group, "dset3");

    // no subgroups: begin and end iterators over subgroups should be equal
    BOOST_TEST_MESSAGE("  empty subgroup container");
    BOOST_CHECK_NO_THROW(groups.begin() == groups.end());

    // use iterator to walk over datasets
    BOOST_TEST_MESSAGE("  walk over datasets");

    container<dataset>::iterator dset_iter = datasets.begin();
    BOOST_CHECK(dset_iter != datasets.end());               // not empty
    BOOST_CHECK(dset_iter->valid());
    BOOST_CHECK_EQUAL(dset_iter.get_name(), "dset1");       // name of element
    BOOST_CHECK_EQUAL(get_name(*dset_iter), "/dset1");      // absolute path

    BOOST_TEST_MESSAGE("  pre-increment");
    BOOST_CHECK(++dset_iter != datasets.end());             // more than 1 element in group
    BOOST_CHECK_EQUAL(dset_iter.get_name(), "dset2");

    BOOST_TEST_MESSAGE("  post-increment");
    BOOST_CHECK_NO_THROW(dset_iter++);
    BOOST_CHECK(dset_iter->valid());
    BOOST_CHECK_EQUAL(dset_iter.get_name(), "dset3");
    BOOST_CHECK(++dset_iter == datasets.end());
}

BOOST_AUTO_TEST_CASE( only_subgroups )
{
    BOOST_TEST_MESSAGE("\nTesting iterator over group that has only subgroups");
    group container_group(file);

    group(container_group, "grp1");
    group(container_group, "grp2");

    container<dataset> datasets = container_group.datasets();
    container<group> groups = container_group.groups();

    // add one more after constructing the container adapter
    group(container_group, "grp3");

    // no datasets: begin and end iterators over datasets should be equal
    BOOST_TEST_MESSAGE("  empty dataset container");
    BOOST_CHECK_NO_THROW(datasets.begin() == datasets.end());

    // use iterator to walk over subgroups
    BOOST_TEST_MESSAGE("  walk over subgroups");

    container<group>::iterator sgroup_iter = groups.begin();
    BOOST_CHECK(sgroup_iter != groups.end());              // not empty
    BOOST_CHECK(sgroup_iter->valid());
    BOOST_CHECK_EQUAL(sgroup_iter.get_name(), "grp1");     // name of element
    BOOST_CHECK_EQUAL(get_name(*sgroup_iter), "/grp1");    // absolute path

    BOOST_TEST_MESSAGE("  post-increment");
    BOOST_CHECK_NO_THROW(sgroup_iter++);
    BOOST_CHECK(sgroup_iter != groups.end());              // more than 1 element in group
    BOOST_CHECK(sgroup_iter->valid());
    BOOST_CHECK_EQUAL(sgroup_iter.get_name(), "grp2");

    BOOST_TEST_MESSAGE("  pre-increment");
    BOOST_CHECK_EQUAL((++sgroup_iter).get_name(), "grp3");
    BOOST_CHECK(++sgroup_iter == groups.end());
}

BOOST_AUTO_TEST_CASE( mixed_group )
{
    BOOST_TEST_MESSAGE("\nTesting iterator over group that both datasets and subgroups");
    group container_group(file);

    create_dataset<int>(container_group, "dset1");
    create_dataset<int>(container_group, "dset2");
    group(container_group, "grp1");
    group(container_group, "grp2");
    group(container_group, "grp3");

    container<dataset> datasets = container_group.datasets();
    container<group> groups = container_group.groups();

    // use dataset iterator to walk over datasets
    BOOST_TEST_MESSAGE("  walk over datasets using post-increment");

    container<dataset>::iterator dset_iter = datasets.begin();
    BOOST_CHECK(dset_iter != datasets.end());               // not empty
    BOOST_CHECK_EQUAL((dset_iter++).get_name(), "dset1");
    BOOST_CHECK_EQUAL((dset_iter++).get_name(), "dset2");
    BOOST_CHECK(dset_iter == datasets.end());

    // use iterator to walk over subgroups
    BOOST_TEST_MESSAGE("  walk over subgroups using pre-increment");

    container<group>::iterator sgroup_iter = groups.begin();
    BOOST_CHECK(sgroup_iter != groups.end());              // not empty
    BOOST_CHECK_EQUAL(sgroup_iter.get_name(), "grp1");
    BOOST_CHECK_EQUAL((++sgroup_iter).get_name(), "grp2");
    BOOST_CHECK_EQUAL((++sgroup_iter).get_name(), "grp3");
    BOOST_CHECK(++sgroup_iter == groups.end());
}

} // namespace fixture
