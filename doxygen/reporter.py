#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

import sys

"""Central issue reporting implementation.

This module implements a Reporter class that is used by other Python modules in
this directory to report issues.  This allows central customization of the
output format, and also a central implementation for redirecting/copying
the output into a log file.  This class also implements sorting for the
messages such that all issues from a single file are reported next to each
other in the output.
"""

class Message(object):

    """Single reported message.

    This class stores the contents of a reporter message for later output to
    allow sorting the output messages reasonably by the reported location.
    """

    def __init__(self, message, details=None, filename=None, location=None):
        """Create a message object.

        The message parameter provides the actual text, while optional details
        provides a list of extra lines that provide context information for the
        error.  filename and location provide two alternative ways of
        specifying the location of the issue:
         - if filename is provided, the issue is reported in that file, without
           a line number
         - if location is provided, get_path() and get_line() methods are
           called to get the file name and line for the message
        """
        if filename:
            self.filename = filename
            self.line = None
        elif location:
            self.filename = location.get_path()
            self.line = location.get_line()
        else:
            self.filename = None
            self.line = None
        self.message = message
        self.details = details

    def __cmp__(self, other):
        """Sort messages based on file name and line number."""
        result = cmp(self.filename, other.filename)
        if not self.filename or result != 0:
            return result
        return cmp(self.line, other.line)

class Reporter(object):

    """Collect and write out issues found by checker scripts."""

    def __init__(self, logfile=None):
        """Initialize the reporter.

        If logfile is set to a file name, all issues will be written to this
        file in addition to stderr.
        """
        self._logfp = None
        if logfile:
            self._logfp = open(logfile, 'w')
        self._messages = []

    def _write(self, message):
        """Implement actual message writing."""
        wholemsg = ''
        if message.filename:
            wholemsg += message.filename + ':'
            if message.line:
                wholemsg += str(message.line) + ':'
            wholemsg += ' '
        wholemsg += message.message
        if message.details:
            wholemsg += '\n    ' + '\n    '.join(message.details)
        wholemsg += '\n'
        sys.stderr.write(wholemsg)
        if self._logfp:
            self._logfp.write(wholemsg)

    def _report(self, message):
        """Handle a single reporter message."""
        if not message.filename:
            self._write(message)
        else:
            self._messages.append(message)

    def write_pending(self):
        """Write out pending messages in sorted order."""
        self._messages.sort()
        for message in self._messages:
            self._write(message)
        self._messages = []

    def close_log(self):
        """Close the log file if one exists."""
        assert not self._messages
        if self._logfp:
            self._logfp.close()
            self._logfp = None

    def xml_assert(self, xmlpath, message):
        """Report issues in Doxygen XML that violate assumptions the script has."""
        self._report(Message('warning: ' + message, filename=xmlpath))

    def input_error(self, message):
        """Report issues in input files."""
        self._report(Message('error: ' + message))

    def file_error(self, fileobj, message):
        """Report file-level issues."""
        self._report(Message('error: ' + message, filename=fileobj.get_path()))

    def code_issue(self, entity, message, details=None):
        """Report an issue in a code construct (not documentation related)."""
        self._report(Message('warning: ' + message, details, location=entity.get_location()))

    def doc_error(self, entity, message):
        """Report an issue in documentation."""
        self._report(Message('error: ' + entity.get_name() + ': ' + message, location=entity.get_location()))

    def doc_note(self, entity, message):
        """Report a potential issue in documentation."""
        self._report(Message('note: ' + entity.get_name() + ': ' + message, location=entity.get_location()))
