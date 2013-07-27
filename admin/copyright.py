#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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
parser.add_option('--check', action='store_true',
                  help='Do not modify the files, only check the copyright. ' +
                       'If specified together with --update, do the modifications ' +
                       'but produce output as if only --check was provided.')
parser.add_option('--update-year', action='store_true',
                  help='Update the copyright year if outdated')
parser.add_option('--update-header', action='store_true',
                  help='Update the copyright header if outdated')
parser.add_option('--add-missing', action='store_true',
                  help='Add missing copyright headers')
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
if not options.check and not options.update_year and \
        not options.update_header and not options.add_missing:
    options.check = True

# The copyright notice to add if it is not present.
class CopyrightNotice(object):
    header = ["This file is part of the GROMACS molecular simulation package.", ""]
    copyright = "Copyright (c) {0} by the GROMACS development team, led by"
    footer = """
Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
and including many others, as listed in the AUTHORS file in the
top-level source directory and at http://www.gromacs.org.

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
""".strip().splitlines()

def guess_filetype(filename):
    basename = os.path.basename(filename)
    root, ext = os.path.splitext(basename)
    if ext == '.cmakein':
        dummy, ext2 = os.path.splitext(root)
        if ext2:
            ext = ext2
    if ext in ('.c', '.cpp', '.h', '.y', '.l'):
        return 'c'
    elif basename == 'CMakeLists.txt' or \
            ext in ('.cmake', '.cmakein', '.py', '.sh'):
        return 'sh'
    return None

def extract_first_comment_block(content_lines, filetype, comment_start):
    if not content_lines or content_lines[0].rstrip() != comment_start:
        return ([], 0)
    comment_block = []
    line_index = 1
    if filetype == 'c':
        while line_index < len(content_lines) and '*/' not in content_lines[line_index]:
            comment_block.append(content_lines[line_index].lstrip('* ').rstrip())
            line_index += 1
    elif filetype == 'sh':
        while line_index < len(content_lines) and content_lines[line_index].startswith('#'):
            comment_block.append(content_lines[line_index].lstrip('# ').rstrip())
            line_index += 1
    return (comment_block, line_index + 1)

def check_copyright(comment_block, copyright_notice):
    copyright_re = r'Copyright \(c\) (([0-9]{4},)+) by the GROMACS development team'
    is_copyright = False
    is_correct = True
    next_header_line = 0
    next_footer_line = 0
    existing_years = ''
    other_copyrights = []
    for line in comment_block:
        if line.startswith('Copyright'):
            is_copyright = True
            match = re.match(copyright_re, line)
            if match:
                existing_years = match.group(1)
            else:
                other_copyrights.append(line)
            if next_header_line != -1:
                is_correct = False
            continue
        if next_header_line >= 0:
            if line == copyright_notice.header[next_header_line]:
                next_header_line += 1
                if next_header_line >= len(copyright_notice.header):
                    next_header_line = -1
            else:
                is_correct = False
        elif next_footer_line >= 0:
            if line == copyright_notice.footer[next_footer_line]:
                next_footer_line += 1
                if next_footer_line >= len(copyright_notice.footer):
                    next_footer_line = -1
            else:
                is_correct = False
        else:
            is_correct = False
    if next_header_line != -1 or next_footer_line != -1:
        is_correct = False

    return (is_copyright, is_correct, existing_years, other_copyrights)

copyright_notice = CopyrightNotice()

# Process each input file in turn.
for filename in filenames:
    # Find the file type from the file name if not explicitly provided.
    filetype = options.lang
    if filename != '-' and not filetype:
        filetype = guess_filetype(filename)
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
        elif filename != '-':
            sys.stderr.write("Unsupported input format: {0}\n".format(filename))
        else:
            sys.stderr.write("No file name or file type provided.\n")
        sys.exit(1)

    # Read the input file.  We are doing an in-place operation, so can't
    # operate in pass-through mode.
    inputfile = sys.stdin
    if filename != '-':
        inputfile = open(filename, 'r')
    contents = inputfile.read().splitlines()
    if filename != '-':
        inputfile.close()

    output = []
    # Keep potential #! line and skip it in the check.
    if contents and contents[0].startswith("#!/"):
        output.append(contents[0])
        contents = contents[1:]

    # Analyze the first comment block in the file.
    comment_block, line_count = extract_first_comment_block(contents, filetype, comment_start)
    is_copyright, is_correct, existing_years, other_copyrights = check_copyright(comment_block, copyright_notice)
    need_update = False

    if existing_years:
        file_years = existing_years
        if not file_years.endswith(years):
            if options.update_year:
                need_update = True
                file_years += years
            if options.check or not options.update_year:
                sys.stderr.write(filename + ': copyright year outdated\n')
            else:
                sys.stderr.write(filename + ': copyright year added\n')
    else:
        file_years = years

    if not is_copyright:
        if options.add_missing:
            need_update = True
        if options.check or not options.add_missing:
            sys.stderr.write(filename + ': copyright header missing\n')
        elif options.add_missing:
            sys.stderr.write(filename + ': copyright header added\n')
    else:
        if not is_correct:
            if options.update_header:
                need_update = True
            if options.check or not need_update:
                sys.stderr.write(filename + ': copyright header outdated\n')
            else:
                sys.stderr.write(filename + ': copyright header updated\n')

    if need_update:
        # Remove the original comment if it was a copyright comment.
        if is_copyright:
            contents = contents[line_count:]
        # Add the new copyright header, potentially with some information
        # extracted from the original comment block.
        output.append(comment_start)
        for cpline in copyright_notice.header:
            line = comment_line + ' ' + cpline
            output.append(line.rstrip())
        for cpline in other_copyrights:
            line = comment_line + ' ' + cpline
            output.append(line.rstrip())
        line = comment_line + ' ' + copyright_notice.copyright.format(file_years)
        output.append(line.rstrip())
        for cpline in copyright_notice.footer:
            line = comment_line + ' ' + cpline
            output.append(line.rstrip())
        output.append(comment_end)

    # Write the output file if required.
    if need_update or filename == '-':
        # Append the rest of the input file as it was.
        output.extend(contents)
        if filename == '-':
            outputfile = sys.stdout
        else:
            outputfile = open(filename, 'w')
        outputfile.write('\n'.join(output) + '\n')
        if filename != '-':
            outputfile.close()
