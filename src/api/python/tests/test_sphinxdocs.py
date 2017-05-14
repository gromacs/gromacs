"""Test any example code embedded in the static module documentation.

I'm not sure if this is completely redundant with sphinx.ext.doctest, but
it is still useful since it will be run with a ``test`` build target and does
not require installation of sphinx. Also, sphinx doctest builds will probably
be separate documentation build target.

Note that ``sphinx-build -b doctest ...`` will use sphinx directives to
configure the test set-up and tear-down, whereas the unit tests will use
the ``setUp`` and ``tearDown`` kwargs to DocFileSuite. For this reason, it
may be best to avoid the need for set-up and tear-down and allow more verbose
examples.

I'm not sure if there is a good way to set more set-up and tear-down from
this script.
"""

import os
import unittest
import doctest

import gmx

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'gmx', 'data'))
gmx.fileio._current_dir = data_dir

def load_tests(loader, tests, ignore):
    # Find sphinx documentation files
    rst_paths = []
    print(data_dir)
    for root, dirs, files in os.walk('../docs'):
        base = os.path.abspath(root)
        rst_paths.extend([os.path.join(base, f) for f in files if f.endswith('.rst')])
    tests.addTests(doctest.DocFileSuite(*rst_paths, module_relative=False))
    return tests
