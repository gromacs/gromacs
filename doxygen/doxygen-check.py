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
from optparse import OptionParser

from doxygengmx import SourceTree, DocType

class Reporter(object):
    def __init__(self, logfile):
        self._logfp = None
        if logfile:
            self._logfp = open(logfile, 'w')

    def _report(self, message):
        sys.stderr.write(message + '\n')
        if self._logfp:
            self._logfp.write(message + '\n')

    def close_log(self):
        if self._logfp:
            self._logfp.close()

    def xml_assert(self, xmlpath, message):
        self._report('warning: {0}: {1}'.format(xmlpath, message))

    def input_error(self, message):
        self._report('error: {0}'.format(message))

    def file_error(self, fileobj, message):
        self._report('error: {0}: {1}'.format(fileobj.get_path(), message))

    def doc_error(self, entity, message):
        self._report('error: {0}: {1}: {2}'.format(entity.get_location(), entity.get_name(), message))

    def doc_note(self, entity, message):
        self._report('note: {0}: {1}: {2}'.format(entity.get_location(), entity.get_name(), message))

def check_file(fileobj, reporter):
    if not fileobj.is_documented():
        # TODO: Add rules for required documentation
        return

    if fileobj.is_source_file():
        # TODO: Add rule to exclude examples from this check
        if fileobj.is_installed():
            reporter.file_error(fileobj, "source file is installed")
        if fileobj.get_documentation_type() != DocType.internal:
            reporter.file_error(fileobj, "source file documentation appears outside full documentation")
        elif fileobj.get_api_type() != DocType.internal:
            reporter.file_error(fileobj, "source file marked as part of an API")
    elif fileobj.is_test_file() and fileobj.is_installed():
        reporter.file_error(fileobj, "test file is installed")
    elif fileobj.is_installed():
        if fileobj.get_documentation_type() != DocType.public:
            reporter.file_error(fileobj, "public header has non-public documentation")
    elif fileobj.get_documentation_type() == DocType.public:
        reporter.file_error(fileobj, "non-installed header has public documentation")
    elif fileobj.get_api_type() == DocType.public:
        reporter.file_error(fileobj, "non-installed header speficied as part of public API")

    if not fileobj.has_brief_description():
        reporter.file_error(fileobj, "is documented, but does not have brief description")

    expectedmod = fileobj.get_expected_module()
    if expectedmod:
        docmodules = fileobj.get_doc_modules()
        if docmodules:
            for module in docmodules:
                if module != expectedmod:
                    reporter.file_error(fileobj, "is documented in incorrect module: {0}".format(module.get_name()))
        else:
            reporter.file_error(fileobj, "is not documented in any module, but {0} exists".format(expectedmod.get_name()))

def check_entity(entity, reporter):
    if entity.is_documented():
        if not entity.has_brief_description():
            reporter.doc_error(entity, "is documented, but does not have brief description")
        if entity.has_inbody_description():
            reporter.doc_note(entity, "has in-body comments, which are ignored")

def main():
    parser = OptionParser()
    parser.add_option('-S', '--source-root',
                      help='Source tree root directory')
    parser.add_option('-B', '--build-root',
                      help='Build tree root directory')
    parser.add_option('--installed',
                      help='Read list of installed files from given file')
    parser.add_option('-l', '--log',
                      help='Write issues into a given log file in addition to stderr')
    options, args = parser.parse_args()

    installedlist = []
    if options.installed:
        with open(options.installed, 'r') as outfile:
            for line in outfile:
                installedlist.append(line.strip())

    reporter = Reporter(options.log)

    tree = SourceTree(options.source_root, options.build_root, reporter)
    tree.set_installed_header_list(installedlist, reporter)

    for fileobj in tree.get_files():
        check_file(fileobj, reporter)

    for entity in tree.get_code_entities():
        check_entity(entity, reporter)

    # TODO: Check #include statements, like old 'make doccheck'

    reporter.close_log()

main()
