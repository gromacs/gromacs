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

class CopyrightState(object):

    """Information about an existing (or non-existing) copyright header."""

    def __init__(self, has_copyright, is_correct, is_newstyle, years, other_copyrights):
        self.has_copyright = has_copyright
        self.is_correct = is_correct
        self.is_newstyle = is_newstyle
        self.years = years
        self.other_copyrights = other_copyrights

class CopyrightChecker(object):

    """Logic for analyzing existing copyright headers and generating new ones."""

    _header = ["This file is part of the GROMACS molecular simulation package.", ""]
    _copyright = "Copyright (c) {0} by the GROMACS development team, led by"
    _footer = """
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

    def check_copyright(self, comment_block):
        """Analyze existing copyright header for correctness and extract information."""
        copyright_re = r'Copyright \(c\) (([0-9]{4},)+) by the GROMACS development team'
        has_copyright = False
        is_newstyle = True
        is_correct = True
        next_header_line = 0
        next_footer_line = 0
        existing_years = ''
        other_copyrights = []
        for line in comment_block:
            if line.startswith('Copyright'):
                has_copyright = True
                match = re.match(copyright_re, line)
                if match:
                    existing_years = match.group(1)
                else:
                    other_copyrights.append(line)
                if next_header_line != -1:
                    is_correct = False
                continue
            if next_header_line >= 0:
                if line == self._header[next_header_line]:
                    next_header_line += 1
                    if next_header_line >= len(self._header):
                        next_header_line = -1
                else:
                    is_correct = False
                    is_newstyle = False
            elif next_footer_line >= 0:
                if line == self._footer[next_footer_line]:
                    next_footer_line += 1
                    if next_footer_line >= len(self._footer):
                        next_footer_line = -1
                else:
                    is_correct = False
            else:
                is_correct = False
        if next_header_line != -1 or next_footer_line != -1:
            is_correct = False

        return CopyrightState(has_copyright, is_correct, is_newstyle, existing_years, other_copyrights)

    def process_copyright(self, state, options, current_years, reporter):
        """Determine whether a copyrigth header needs to be updated and report issues."""
        need_update = False

        if state.years:
            new_years = state.years
            if not new_years.endswith(current_years):
                if options.update_year:
                    need_update = True
                    new_years += current_years
                if options.check or not need_update:
                    reporter.report('copyright year outdated')
                else:
                    reporter.report('copyright year added')
        else:
            new_years = current_years

        if not state.has_copyright:
            if options.add_missing:
                need_update = True
            if options.check or not need_update:
                reporter.report('copyright header missing')
            elif options.add_missing:
                reporter.report('copyright header added')
        else:
            if not state.is_newstyle:
                if options.replace_header:
                    need_update = True
                if options.check or not need_update:
                    reporter.report('copyright header incorrect')
                else:
                    reporter.report('copyright header replaced')
            elif not state.is_correct:
                if options.update_header:
                    need_update = True
                if options.check or not need_update:
                    reporter.report('copyright header outdated')
                else:
                    reporter.report('copyright header updated')

        return need_update, new_years

    def get_copyright_text(self, years, other_copyrights):
        """Construct a new copyright header."""
        output = []
        output.extend(self._header)
        if other_copyrights:
            output.extend(other_copyrights)
        output.append(self._copyright.format(years))
        output.extend(self._footer)
        return output

class Reporter(object):

    """Wrapper for reporting issues in a file."""

    def __init__(self, reportfile, filename):
        self._reportfile = reportfile
        self._filename = filename

    def report(self, text):
        self._reportfile.write(self._filename + ': ' + text + '\n');

class CommentHandlerC(object):

    """Handler for extracting and creating C-style comments."""

    def extract_first_comment_block(self, content_lines):
        if not content_lines or content_lines[0].rstrip() != '/*':
            return ([], 0)
        comment_block = []
        line_index = 1
        while line_index < len(content_lines) and '*/' not in content_lines[line_index]:
            comment_block.append(content_lines[line_index].lstrip('* ').rstrip())
            line_index += 1
        return (comment_block, line_index + 1)

    def create_comment_block(self, lines):
        output = []
        output.append('/*')
        output.extend([(' * ' + x).rstrip() for x in lines])
        output.append(' */')
        return output

class CommentHandlerSh(object):

    """Handler for extracting and creating sh-style comments."""

    def extract_first_comment_block(self, content_lines):
        if not content_lines or content_lines[0].rstrip() != '#':
            return ([], 0)
        comment_block = []
        line_index = 1
        while line_index < len(content_lines) and content_lines[line_index].startswith('#'):
            comment_block.append(content_lines[line_index].lstrip('* ').rstrip())
            line_index += 1
        return (comment_block, line_index + 1)

    def create_comment_block(self, lines):
        output = []
        output.append('#')
        output.extend([('# ' + x).rstrip() for x in lines])
        output.append('')
        return output

comment_handlers = {'c': CommentHandlerC(), 'sh': CommentHandlerSh()}

def select_comment_handler(override, filename):
    """Select comment handler for a file based on file name and input options."""
    filetype = override
    if not filetype and filename != '-':
        basename = os.path.basename(filename)
        root, ext = os.path.splitext(basename)
        if ext == '.cmakein':
            dummy, ext2 = os.path.splitext(root)
            if ext2:
                ext = ext2
        if ext in ('.c', '.cpp', '.h', '.y', '.l', '.pre'):
            filetype = 'c'
        elif basename in ('CMakeLists.txt', 'GMXRC', 'git-pre-commit') or \
                ext in ('.cmake', '.cmakein', '.py', '.sh', '.bash', '.csh', '.zsh'):
            filetype = 'sh'
    if filetype in comment_handlers:
        return comment_handlers[filetype]
    if filetype:
        sys.stderr.write("Unsupported input format: {0}\n".format(filetype))
    elif filename != '-':
        sys.stderr.write("Unsupported input format: {0}\n".format(filename))
    else:
        sys.stderr.write("No file name or file type provided.\n")
    sys.exit(1)

def process_options():
    """Process input options."""
    parser = OptionParser()
    parser.add_option('-l', '--lang',
                      help='Comment type to use (c or sh)')
    parser.add_option('-y', '--years',
                      help='Comma-separated list of years')
    parser.add_option('-F', '--files',
                      help='File to read list of files from')
    parser.add_option('--check', action='store_true',
                      help='Do not modify the files, only check the copyright (default action). ' +
                           'If specified together with --update, do the modifications ' +
                           'but produce output as if only --check was provided.')
    parser.add_option('--update-year', action='store_true',
                      help='Update the copyright year if outdated')
    parser.add_option('--update-header', action='store_true',
                      help='Update the copyright header if outdated')
    parser.add_option('--replace-header', action='store_true',
                      help='Replace any copyright header with the current one')
    parser.add_option('--add-missing', action='store_true',
                      help='Add missing copyright headers')
    options, args = parser.parse_args()

    filenames = args
    if options.files:
        with open(options.files, 'r') as filelist:
            filenames = [x.strip() for x in filelist.read().splitlines()]
    elif not filenames:
        filenames = ['-']

    # Default is --check if nothing provided.
    if not options.check and not options.update_year and \
            not options.update_header and not options.replace_header and \
            not options.add_missing:
        options.check = True

    return options, filenames

# Main execution begins here.

options, filenames = process_options()
years = options.years
if not years:
    years = str(datetime.date.today().year)
if not years.endswith(','):
    years += ','

checker = CopyrightChecker()

# Process each input file in turn.
for filename in filenames:
    comment_handler = select_comment_handler(options.lang, filename)

    # Read the input file.  We are doing an in-place operation, so can't
    # operate in pass-through mode.
    if filename == '-':
        contents = sys.stdin.read().splitlines()
        reporter = Reporter(sys.stderr, '<stdin>')
    else:
        with open(filename, 'r') as inputfile:
            contents = inputfile.read().splitlines()
        reporter = Reporter(sys.stdout, filename)

    output = []
    # Keep potential #! line and skip it in the check.
    if contents and contents[0].startswith("#!/"):
        output.append(contents[0])
        contents = contents[1:]

    # Analyze the first comment block in the file.
    comment_block, line_count = comment_handler.extract_first_comment_block(contents)
    state = checker.check_copyright(comment_block)
    need_update, file_years = checker.process_copyright(state, options, years, reporter)

    if need_update:
        # Remove the original comment if it was a copyright comment.
        if state.has_copyright:
            contents = contents[line_count:]
        new_block = checker.get_copyright_text(file_years, state.other_copyrights)
        output.extend(comment_handler.create_comment_block(new_block))

    # Write the output file if required.
    if need_update or filename == '-':
        # Append the rest of the input file as it was.
        output.extend(contents)
        output = '\n'.join(output) + '\n'
        if filename == '-':
            sys.stdout.write(output)
        else:
            with open(filename, 'w') as outputfile:
                outputfile.write(output)
