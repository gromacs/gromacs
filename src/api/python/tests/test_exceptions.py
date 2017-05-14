"""Test gmx module exceptions

Regression tests for defined exceptions and inheritance.
Throwing of appropriate exceptions is tested using assertRaises
in test for the components that throw them.
"""

import unittest

import gmx
from gmx.exceptions import Error
from gmx.exceptions import UsageError

# Note: should note API level...

class ExceptionTestCase(unittest.TestCase):
    def test_exception_inheritance(self):
        exception = None
        try:
            raise UsageError("generic usage error")
        except gmx.exceptions.UsageError as e:
            exception = e
        self.assertTrue(isinstance(exception, UsageError))
        self.assertTrue(isinstance(exception, Error))
