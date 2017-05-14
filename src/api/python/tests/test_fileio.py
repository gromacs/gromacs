"""Test gmx.fileio submodule"""

import unittest
import os

from gmx.fileio import TprFile
from gmx.exceptions import UsageError

tpr_filename = os.path.join(os.path.dirname(__file__), 'data', 'test.tpr')

class TprTestCase(unittest.TestCase):
    def test_tprfile_read(self):
        self.assertRaises(UsageError, TprFile)
        self.assertRaises(UsageError, TprFile, tpr_filename, 'x')
        # TprFile does not yet check whether file exists and is readable...
        #self.assertRaises(UsageError, TprFile, 1, 'r')
        fh = TprFile(tpr_filename, 'r')
