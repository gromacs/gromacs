"""Test embedded documentaion-by-example in docstrings."""

import unittest
import doctest
import gmx

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(gmx))
    tests.addTests(doctest.DocTestSuite(gmx.exceptions))
    tests.addTests(doctest.DocTestSuite(gmx.fileio))
    tests.addTests(doctest.DocTestSuite(gmx.util))
    return tests
