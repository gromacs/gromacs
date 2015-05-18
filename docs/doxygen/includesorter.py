#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

"""Include directive sorter for GROMACS.

This module implements an #include directive sorter for GROMACS C/C++ files.
It allows (in most cases) automatically sorting includes and formatting
the paths to use either relative paths or paths relative to src/.
It groups includes in groups of related headers, sorts the headers
alphabetically within each block, and inserts empty lines in between.
It can be run as a standalone script, in which case it requires an up-to-date
list of installed headers and Doxygen XML documentation to be present in the
build tree.  It can also be imported as a module to be embedded in other
scripts.  In the latter case, the IncludeSorter provides the main interface.

The sorting assumes some conventions (e.g., that system headers are included
with angle brackets instead of quotes).  Generally, these conventions are
checked by the check-source.py script.
"""

import os.path
import re
import sys

class IncludeGroup(object):

    """Enumeration type for grouping includes."""

    def __init__(self, value):
        """Initialize a IncludeGroup instance.

        IncludeGroup.{main,system_c,...} should be used outside the
        class instead of calling the constructor.
        """
        self._value = value

    def __cmp__(self, other):
        """Order include groups in the desired order."""
        return cmp(self._value, other._value)

# gmxpre.h is always first
IncludeGroup.pre = IncludeGroup(0)
# "main" include file for the source file is next
IncludeGroup.main = IncludeGroup(1)
# config.h is next, if present, to keep its location consistent
IncludeGroup.config = IncludeGroup(2)
# Followed by system headers, with C first and C++ following
IncludeGroup.system_c = IncludeGroup(3)
IncludeGroup.system_c_cpp = IncludeGroup(4)
IncludeGroup.system_cpp = IncludeGroup(5)
# System headers not in standard C/C++ are in a separate block
IncludeGroup.system_other = IncludeGroup(6)
# src/external/ contents that are included with quotes go here
IncludeGroup.nonsystem_other = IncludeGroup(7)
# Other GROMACS headers
IncludeGroup.gmx_general = IncludeGroup(8)
# This group is for shared (unit) testing utilities
IncludeGroup.gmx_test = IncludeGroup(9)
# This group is for headers local to the including file/module
IncludeGroup.gmx_local = IncludeGroup(10)

class GroupedSorter(object):

    """Grouping and formatting logic for #include directives.

    This class implements the actual logic that decides how includes are
    grouped and sorted, and how they are formatted."""

    # These variables contain the list of system headers for various blocks
    _std_c_headers = ['assert.h', 'ctype.h', 'errno.h', 'float.h',
            'inttypes.h', 'limits.h', 'math.h', 'signal.h', 'stdarg.h',
            'stddef.h', 'stdint.h', 'stdio.h', 'stdlib.h', 'string.h',
            'time.h']
    _std_c_cpp_headers = ['c' + x[:-2] for x in _std_c_headers]
    _std_cpp_headers = ['algorithm', 'deque', 'exception', 'fstream',
            'iomanip', 'ios', 'iosfwd', 'iostream', 'istream', 'iterator',
            'limits', 'list', 'map', 'memory', 'new', 'numeric', 'ostream',
            'regex', 'set', 'sstream', 'stdexcept', 'streambuf', 'string', 'strstream',
            'typeinfo', 'vector', 'utility']

    def __init__(self, style='pub-priv', absolute=False):
        """Initialize a sorted with the given style."""
        if style == 'single-group':
            self._local_group = 'none'
        elif style == 'pub-priv':
            self._local_group = 'private'
        else:
            self._local_group = 'local'
        if absolute:
            self._abspath_main = True
            self._abspath_local = True
        else:
            self._abspath_main = False
            self._abspath_local = False

    def _get_path(self, included_file, group, including_file):
        """Compute include path to use for an #include.

        The path is made either absolute (i.e., relative to src/), or
        relative to the location of the including file, depending on the group
        the file is in.
        """
        use_abspath = including_file is None or group is None
        if not use_abspath:
            if group in (IncludeGroup.gmx_general, IncludeGroup.gmx_test):
                use_abspath = True
            elif group == IncludeGroup.main and self._abspath_main:
                use_abspath = True
            elif group == IncludeGroup.gmx_local and self._abspath_local:
                use_abspath = True
        if not use_abspath:
            fromdir = os.path.dirname(including_file.get_abspath())
            relpath = os.path.relpath(included_file.get_abspath(), fromdir)
            if not relpath.startswith('..'):
                return relpath
        path = included_file.get_relpath()
        assert path.startswith('src/')
        return path[4:]

    def _get_gmx_group(self, including_file, included_file):
        """Determine group for GROMACS headers.

        Helper function to determine the group for an #include directive
        when the #include is in one of the gmx_* groups (or in the main group).
        """
        main_header = including_file.get_main_header()
        if main_header and main_header == included_file:
            return IncludeGroup.main
        if included_file.get_directory().get_name() == 'testutils':
            return IncludeGroup.gmx_test
        if including_file.get_directory().contains(included_file):
            if self._local_group == 'local':
                return IncludeGroup.gmx_local
            if self._local_group == 'private':
                if included_file.api_type_is_reliable() \
                        and included_file.is_module_internal():
                    return IncludeGroup.gmx_local
                if not included_file.api_type_is_reliable() \
                        and including_file.get_relpath().startswith('src/programs'):
                    return IncludeGroup.gmx_local
        if included_file.is_test_file():
            return IncludeGroup.gmx_test
        return IncludeGroup.gmx_general

    def get_sortable_object(self, include):
        """Produce a sortable, opaque object for an include.

        Includes are sorted by calling this function for each #include object,
        and sorting the list made up of these objects (using the default
        comparison operators).  Each element from the sorted list is then
        passed to format_include(), which extracts information from the opaque
        object and formats the #include directive for output.
        """
        included_file = include.get_file()
        if not included_file:
            path = include.get_included_path()
            if path in self._std_c_headers:
                group = IncludeGroup.system_c
            elif path in self._std_c_cpp_headers:
                group = IncludeGroup.system_c_cpp
            elif path in self._std_cpp_headers:
                group = IncludeGroup.system_cpp
            else:
                group = IncludeGroup.system_other
        elif included_file.is_external():
            group = IncludeGroup.nonsystem_other
            if 'external/' in include.get_included_path():
                path = self._get_path(included_file, group, None)
            else:
                path = include.get_included_path()
        elif included_file.get_name() == 'gmxpre.h':
            group = IncludeGroup.pre
            path = self._get_path(included_file, group, None)
        elif included_file.get_name() == 'config.h':
            group = IncludeGroup.config
            path = self._get_path(included_file, group, None)
        else:
            including_file = include.get_including_file()
            group = self._get_gmx_group(including_file, included_file)
            path = self._get_path(included_file, group, including_file)
        return (group, os.path.split(path), include)

    def format_include(self, obj, prev):
        """Format an #include directive after sorting."""
        result = []
        if prev:
            if prev[0] != obj[0]:
                # Print empty line between groups
                result.append('\n')
            elif prev[1] == obj[1]:
                # Skip duplicates
                return result
        include = obj[2]
        line = include.get_full_line()
        include_re = r'^(?P<head>\s*#\s*include\s+)["<][^">]*[">](?P<tail>.*)$'
        match = re.match(include_re, line)
        assert match
        if include.is_system():
            path = '<{0}>'.format(os.path.join(obj[1][0], obj[1][1]))
        else:
            path = '"{0}"'.format(os.path.join(obj[1][0], obj[1][1]))
        result.append('{0}{1}{2}\n'.format(match.group('head'), path, match.group('tail')))
        return result

class IncludeSorter(object):

    """High-level logic for sorting includes.

    This class contains the high-level logic for sorting include statements.
    The actual ordering and formatting the includes is delegated to a sort method
    (see GroupedSorter) to keep things separated.
    """

    def __init__(self, sortmethod=None, quiet=True):
        """Initialize the include sorter with the given sorter and options."""
        if not sortmethod:
            sortmethod = GroupedSorter()
        self._sortmethod = sortmethod
        self._quiet = quiet
        self._changed = False

    def _sort_include_block(self, block, lines):
        """Sort a single include block.

        Returns a new list of lines for the block.
        If anything is changed, self._changed is set to True, and the caller
        can check that."""
        includes = map(self._sortmethod.get_sortable_object, block.get_includes())
        includes.sort()
        result = []
        prev = None
        current_line_number = block.get_first_line()-1
        for include in includes:
            newlines = self._sortmethod.format_include(include, prev)
            result.extend(newlines)
            if not self._changed:
                for offset, newline in enumerate(newlines):
                    if lines[current_line_number + offset] != newline:
                        self._changed = True
                        break
                current_line_number += len(newlines)
            prev = include
        return result

    def sort_includes(self, fileobj):
        """Sort all includes in a file."""
        lines = fileobj.get_contents()
        # Format into a list first:
        #  - avoid bugs or issues in the script truncating the file
        #  - can check whether anything was changed before touching the file
        newlines = []
        prev = 0
        self._changed = False
        for block in fileobj.get_include_blocks():
            newlines.extend(lines[prev:block.get_first_line()-1])
            newlines.extend(self._sort_include_block(block, lines))
            # The returned values are 1-based, but indexing here is 0-based,
            # so an explicit +1 is not needed.
            prev = block.get_last_line()
        if self._changed:
            if not self._quiet:
                sys.stderr.write('{0}: includes reformatted\n'.format(fileobj.get_relpath()))
            newlines.extend(lines[prev:])
            with open(fileobj.get_abspath(), 'w') as fp:
                fp.write(''.join(newlines))

    def check_sorted(self, fileobj):
        """Check that includes within a file are sorted."""
        # TODO: Make the checking work without full contents of the file
        lines = fileobj.get_contents()
        is_sorted = True
        details = None
        for block in fileobj.get_include_blocks():
            self._changed = False
            sorted_lines = self._sort_include_block(block, lines)
            if self._changed:
                is_sorted = False
                # TODO: Do a proper diff to show the actual changes.
                if details is None:
                    details = ["Correct order/style is:"]
                else:
                    details.append("    ...")
                details.extend(["    " + x.rstrip() for x in sorted_lines])
        return (is_sorted, details)

def main():
    """Run the include sorter script."""
    import os
    import sys

    from optparse import OptionParser

    from gmxtree import GromacsTree
    from reporter import Reporter

    parser = OptionParser()
    parser.add_option('-S', '--source-root',
                      help='Source tree root directory')
    parser.add_option('-B', '--build-root',
                      help='Build tree root directory')
    parser.add_option('-F', '--files',
                      help='Specify files to sort')
    parser.add_option('-q', '--quiet', action='store_true',
                      help='Do not write status messages')
    # This is for evaluating different options; can be removed from the final
    # version.
    parser.add_option('-s', '--style', type='choice', default='pub-priv',
                      choices=('single-group', 'pub-priv', 'pub-local'),
                      help='Style for GROMACS includes')
    parser.add_option('--absolute', action='store_true',
                      help='Write all include paths relative to src/')
    options, args = parser.parse_args()

    filelist = args
    if options.files:
        if options.files == '-':
            lines = sys.stdin.readlines()
        else:
            with open(options.files, 'r') as fp:
                lines = fp.readlines()
        filelist.extend([x.strip() for x in lines])

    reporter = Reporter(quiet=True)

    if not options.quiet:
        sys.stderr.write('Scanning source tree...\n')
    tree = GromacsTree(options.source_root, options.build_root, reporter)
    tree.load_installed_file_list()
    files = []
    for filename in filelist:
        fileobj = tree.get_file(os.path.abspath(filename))
        if not fileobj:
            sys.stderr.write('warning: ignoring unknown file {0}\n'.format(filename))
            continue
        files.append(fileobj)
    if not options.quiet:
        sys.stderr.write('Reading source files...\n')
    tree.scan_files(only_files=files, keep_contents=True)
    extfiles = set(files)
    for fileobj in files:
        for included_file in fileobj.get_includes():
            other_file = included_file.get_file()
            if other_file:
                extfiles.add(other_file)
    if not options.quiet:
        sys.stderr.write('Reading Doxygen XML files...\n')
    tree.load_xml(only_files=extfiles)

    if not options.quiet:
        sys.stderr.write('Sorting includes...\n')

    sorter = IncludeSorter(GroupedSorter(options.style, options.absolute), options.quiet)

    for fileobj in files:
        sorter.sort_includes(fileobj)

if __name__ == '__main__':
    main()
