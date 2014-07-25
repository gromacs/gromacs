#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

import os.path
import re

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

IncludeGroup.pre = IncludeGroup(0)
IncludeGroup.main = IncludeGroup(1)
IncludeGroup.config = IncludeGroup(2)
IncludeGroup.system_c = IncludeGroup(3)
IncludeGroup.system_c_cpp = IncludeGroup(4)
IncludeGroup.system_cpp = IncludeGroup(5)
IncludeGroup.system_other = IncludeGroup(6)
IncludeGroup.nonsystem_other = IncludeGroup(7)
IncludeGroup.gmx_general = IncludeGroup(8)
IncludeGroup.gmx_test = IncludeGroup(9)
IncludeGroup.gmx_local = IncludeGroup(10)

class GroupedSorter(object):
    # TODO: Ensure that these are complete enough.
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

    def __init__(self, style):
        if style == 'single-group':
            self._local_group = 'none'
        elif style.startswith('pub-priv'):
            self._local_group = 'private'
        else:
            self._local_group = 'local'
        if style.endswith('-rel'):
            self._abspath_main = False
            self._abspath_local = False
        else:
            self._abspath_main = True
            self._abspath_local = True

    def _get_path(self, included_file, group, including_file):
        use_abspath = including_file is None or group is None
        if not use_abspath:
            if group == IncludeGroup.gmx_general or group == IncludeGroup.gmx_test:
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
        elif included_file.get_name() in ('config.h', 'gmx_header_config.h'):
            group = IncludeGroup.config
            path = self._get_path(included_file, group, None)
        else:
            including_file = include.get_including_file()
            group = self._get_gmx_group(including_file, included_file)
            path = self._get_path(included_file, group, including_file)
        return (group, os.path.split(path), include)

    def format_include(self, obj, prev, lines):
        result = []
        if prev:
            if prev[0] != obj[0]:
                # Print empty line between groups
                result.append('\n')
            elif prev[1] == obj[1]:
                # Skip duplicates
                return result
        include = obj[2]
        line = lines[include.get_line_number()-1]
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
    def __init__(self, sortmethod):
        self._sortmethod = sortmethod

    def _sort_include_block(self, block, lines):
        includes = map(self._sortmethod.get_sortable_object, block.get_includes())
        includes.sort()
        result = []
        prev = None
        for include in includes:
            result.extend(self._sortmethod.format_include(include, prev, lines))
            prev = include
        return result

    def sort_includes(self, fileobj):
        lines = fileobj.get_contents()
        # Format into a list first to avoid bugs or issues in the script
        # truncating the file.
        newlines = []
        prev = 0
        for block in fileobj.get_include_blocks():
            newlines.extend(lines[prev:block.get_first_line()-1])
            newlines.extend(self._sort_include_block(block, lines))
            # The returned values are 1-based, but indexing here is 0-based,
            # so an explicit +1 is not needed.
            prev = block.get_last_line()
        newlines.extend(lines[prev:])
        with open(fileobj.get_abspath(), 'w') as fp:
            fp.write(''.join(newlines))

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
    parser.add_option('-s', '--style', type='choice', default='pub-priv-rel',
                      choices=('single-group', 'pub-priv-abs', 'pub-priv-rel',
                          'pub-local-abs', 'pub-local-rel'),
                      help='Style for Gromacs includes')
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

    sorter = IncludeSorter(GroupedSorter(options.style))

    for fileobj in files:
        sorter.sort_includes(fileobj)

if __name__ == '__main__':
    main()
