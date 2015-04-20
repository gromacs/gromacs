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

from fnmatch import fnmatch

"""Central issue reporting implementation.

This module implements a Reporter class that is used by other Python modules in
this directory to report issues.  This allows central customization of the
output format, and also a central implementation for redirecting/copying
the output into a log file.  This class also implements sorting for the
messages such that all issues from a single file are reported next to each
other in the output, as well as filtering to make it possible to suppress
certain messages.
"""

class Location(object):

    """Location for a reported message."""

    def __init__(self, filename, line):
        """Create a location with the given file and line number.

        One or both of the parameters can be None, but filename should be
        specified if line is.
        """
        self.filename = filename
        self.line = line

    def __nonzero__(self):
        """Make empty locations False in boolean context."""
        return self.filename is not None

    def __str__(self):
        """Format the location as a string."""
        if self.line:
            return '{0}:{1}'.format(self.filename, self.line)
        elif self.filename:
            return self.filename
        else:
            return '<unknown>'

    def __cmp__(self, other):
        """Sort locations based on file name and line number."""
        result = cmp(self.filename, other.filename)
        if not self.filename or result != 0:
            return result
        return cmp(self.line, other.line)

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
         - if location is provided, it should be a Location instance
        """
        if filename:
            self.location = Location(filename, None)
        elif location:
            self.location = location
        else:
            self.location = Location(None, None)
        self.message = message
        self.details = details

    def __cmp__(self, other):
        """Sort messages based on file name and line number."""
        return cmp(self.location, other.location)

class Filter(object):

    """Filter expression to exclude messages."""

    def __init__(self, filterline):
        """Initialize a filter from a line in a filter file."""
        self._orgline = filterline
        filepattern, text = filterline.split(':', 1)
        if filepattern == '*':
            self._filematcher = lambda x: x is not None
        elif filepattern:
            self._filematcher = lambda x: x and fnmatch(x, '*/' + filepattern)
        else:
            self._filematcher = lambda x: x is None
        self._textpattern = text.strip()
        self._count = 0

    def matches(self, message):
        """Check whether the filter matches a message."""
        if not self._filematcher(message.location.filename):
            return False
        if not fnmatch(message.message, self._textpattern):
            return False
        self._count += 1
        return True

    def get_match_count(self):
        """Return the number of times this filter has matched."""
        return self._count

    def get_text(self):
        """Return original line used to specify the filter."""
        return self._orgline

class Reporter(object):

    """Collect and write out issues found by checker scripts."""

    def __init__(self, logfile=None, quiet=False):
        """Initialize the reporter.

        If logfile is set to a file name, all issues will be written to this
        file in addition to stderr.

        If quiet is set to True, the reporter will suppress all output.
        """
        self._logfp = None
        if logfile:
            self._logfp = open(logfile, 'w')
        self._messages = []
        self._filters = []
        self._quiet = quiet
        self._had_warnings = False

    def _write(self, message):
        """Implement actual message writing."""
        wholemsg = ''
        if message.location:
            wholemsg += str(message.location) + ': '
        wholemsg += message.message
        if message.details:
            wholemsg += '\n    ' + '\n    '.join(message.details)
        wholemsg += '\n'
        sys.stderr.write(wholemsg)
        if self._logfp:
            self._logfp.write(wholemsg)
        self._had_warnings = True

    def _report(self, message):
        """Handle a single reporter message."""
        if self._quiet:
            return
        for filterobj in self._filters:
            if filterobj.matches(message):
                return
        if not message.location:
            self._write(message)
        else:
            self._messages.append(message)

    def load_filters(self, filterfile):
        """Load filters for excluding messages from a file."""
        with open(filterfile, 'r') as fp:
            for filterline in fp:
                filterline = filterline.strip()
                if not filterline or filterline.startswith('#'):
                    continue
                self._filters.append(Filter(filterline))

    def write_pending(self):
        """Write out pending messages in sorted order."""
        self._messages.sort()
        for message in self._messages:
            self._write(message)
        self._messages = []

    def report_unused_filters(self):
        """Report filters that did not match any messages."""
        for filterobj in self._filters:
            if filterobj.get_match_count() == 0:
                # TODO: Consider adding the input filter file as location
                text = 'warning: unused filter: ' + filterobj.get_text()
                self._write(Message(text))

    def had_warnings(self):
        """Return true if any warnings have been reported."""
        return self._had_warnings

    def close_log(self):
        """Close the log file if one exists."""
        assert not self._messages
        if self._logfp:
            self._logfp.close()
            self._logfp = None

    def xml_assert(self, xmlpath, message):
        """Report issues in Doxygen XML that violate assumptions in the script."""
        self._report(Message('warning: ' + message, filename=xmlpath))

    def input_error(self, message):
        """Report issues in input files."""
        self._report(Message('error: ' + message))

    def file_error(self, fileobj, message):
        """Report file-level issues."""
        self._report(Message('error: ' + message,
            location=fileobj.get_reporter_location()))

    def code_issue(self, entity, message, details=None):
        """Report an issue in a code construct (not documentation related)."""
        self._report(Message('warning: ' + message, details,
            location=entity.get_reporter_location()))

    def cyclic_issue(self, message, details=None):
        """Report a cyclic dependency issue."""
        self._report(Message('warning: ' + message, details))

    def doc_error(self, entity, message):
        """Report an issue in documentation."""
        self._report(Message('error: ' + entity.get_name() + ': ' + message,
            location=entity.get_reporter_location()))

    def doc_note(self, entity, message):
        """Report a potential issue in documentation."""
        self._report(Message('note: ' + entity.get_name() + ': ' + message,
            location=entity.get_reporter_location()))
