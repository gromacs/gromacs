#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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

import datetime
import os.path
import re
import sys

from optparse import OptionParser

# Declare options.
parser = OptionParser()
parser.add_option('-l', '--lang',
                  help='Comment type to use (c or sh)')
parser.add_option('-y', '--years',
                  help='Comma-separated list of years')
parser.add_option('-F', '--files',
                  help='File to read list of files from')
parser.add_option('--input-suffix',
                  help='Suffix to add to all provided input files')
parser.add_option('--check', action='store_true',
                  help='Do not modify the files, only check the copyright. ' +
                       'If specified together with --update, do the modifications ' +
                       'but produce output as if only --check was provided.')
parser.add_option('--update', action='store_true',
                  help='Modify the source files (default action)')
parser.add_option('--no-add', action='store_true',
                  help='Do not add missing copyright headers')
options, args = parser.parse_args()

# Process options.
filenames = args
if options.files:
    with open(options.files, 'r') as filelist:
        filenames = [x.strip() for x in filelist.read().splitlines()]
elif not filenames:
    filenames = ['-']
years = options.years
if not years:
    years = str(datetime.date.today().year)
if not years.endswith(','):
    years += ','
update = options.update or not options.check

# The copyright notice to add if it is not present.
copyright_notice_lines = """
This file is part of the GROMACS molecular simulation package.

Copyright (c) {0} by the GROMACS development team, led by
David van der Spoel, Berk Hess, Erik Lindahl, and including many
others, as listed in the AUTHORS file in the top-level source
directory and at http://www.gromacs.org.

GROMACS is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 2.1
of the License, or (at your option) any later version.

GROMACS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with GROMACS; if not, see
http://www.gnu.org/licenses, or write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.

If you want to redistribute modifications to GROMACS, please
consider that scientific software is very special. Version
control is crucial - bugs must be traceable. We will be happy to
consider code for inclusion in the official distribution, but
derived work must not be called official GROMACS. Details are found
in the README & COPYING files - if they are missing, get the
official version at http://www.gromacs.org.

To help us fund GROMACS development, we humbly ask that you cite
the research papers on the package. Check out http://www.gromacs.org.
""".format(years).strip().splitlines()

# Process each input file in turn.
for filename in filenames:
    # Find the file type from the file name if not explicitly provided.
    filetype = options.lang
    if filename != '-' and not filetype:
        basename = os.path.basename(filename)
        root, ext = os.path.splitext(basename)
        if ext == '.cmakein':
            dummy, ext2 = os.path.splitext(root)
            if ext2:
                ext = ext2
        if ext in ('.c', '.cpp', '.h', '.y', '.l'):
            filetype = 'c'
        elif basename == 'CMakeLists.txt' or \
                ext in ('.cmake', '.cmakein', '.py', '.sh'):
            filetype = 'sh'
    # Determine the comment characters from the file type.
    if filetype == 'c':
        comment_start = "/*"
        comment_line = " *"
        comment_end = " */"
    elif filetype == 'sh':
        comment_start = "#"
        comment_line = "#"
        comment_end = ""
    else:
        if filetype:
            sys.stderr.write("Unsupported input format: {0}\n".format(filetype))
        elif filename:
            sys.stderr.write("Unsupported input format: {0}\n".format(filename))
        else:
            sys.stderr.write("No file name or file type provided.\n")
        sys.exit(1)

    # Read the input file.  We may be doing an in-place operation, so can't
    # operate in pass-through mode.
    inputfile = sys.stdin
    if filename != '-':
        if options.input_suffix:
            filename = filename + options.input_suffix
        inputfile = open(filename, 'r')
    contents = inputfile.read().splitlines()
    if filename != '-':
        inputfile.close()

    # Open the output file if required
    if not update:
        outputfile = None
    elif filename == '-':
        outputfile = sys.stdout
    else:
        outputfile = open(filename, 'w')
    # Keep potential #! line and skip it in the check.
    if contents[0].startswith("#!/"):
        if outputfile:
            outputfile.write(contents[0] + '\n')
        contents = contents[1:]
    # Check whether the file starts with the copyright notice
    # (only check the first content line)
    if contents[0].rstrip() == comment_start and \
            contents[1].rstrip() == comment_line + ' ' + copyright_notice_lines[0].rstrip():
        # regexp to find the copyright line where to check/update the years.
        copyright_re = re.escape(comment_line) + r' Copyright \(c\) (([0-9]{4},)+) by the GROMACS development team, led by'
        # Check the copyright notice until the end of the first comment.
        # TODO: Consider whether also the subsequent lines of the notice should
        # be checked for correctness.
        for linenr in range(1, len(contents)):
            line = contents[linenr]
            if not line.startswith(comment_line) or line.rstrip() == comment_end:
                break

            # Check that the year is present on the correct line.
            # TODO: Deal with the issue if no line matches the regexp.
            match = re.match(copyright_re, line)
            if match:
                if not match.group(1).endswith(years):
                    if options.check:
                        sys.stderr.write(filename + ': copyright year outdated\n')
                    else:
                        sys.stderr.write(filename + ': copyright year added\n')
                    if update:
                        pos = match.end(1)
                        contents[linenr] = line[:pos] + years + line[pos:]
    else:
        # The file does not start with the correct header.
        if options.check or options.no_add:
            sys.stderr.write(filename + ': copyright header missing\n')
        else:
            sys.stderr.write(filename + ': copyright header added\n')
        if update and not options.no_add:
            outputfile.write(comment_start + '\n')
            for cpline in copyright_notice_lines:
                line = comment_line + ' ' + cpline
                outputfile.write(line.rstrip() + '\n')
            outputfile.write(comment_end + '\n')
    # Write the rest of the output file as it was.
    if update:
        outputfile.write('\n'.join(contents) + '\n')
        if filename and filename != '-':
            outputfile.close()
