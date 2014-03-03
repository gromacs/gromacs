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

import os.path
import sys

import doxygenxml as xml

class DocType(object):
    none = 'none'
    public = 'public'
    library = 'library'
    internal = 'internal'

def _get_api_type_for_compound(grouplist):
    result = DocType.internal
    for group in grouplist:
        if isinstance(group, xml.Group):
            if group.get_name() == 'group_publicapi':
                result = DocType.public
            elif group.get_name() == 'group_libraryapi':
                result = DocType.library
            # TODO: Check for multiple group membership
    return result

class File(object):
    def __init__(self, rawdoc):
        self._rawdoc = rawdoc
        self._installed = False
        dummy, extension = os.path.splitext(rawdoc.get_path())
        self._sourcefile = False
        if extension in ('.c', '.cpp', '.cu'):
            self._sourcefile = True
        self._dir = None
        self._apitype = DocType.none
        if self._rawdoc.is_documented():
            self._apitype = _get_api_type_for_compound(self._rawdoc.get_groups())

    def set_directory(self, dirobj):
        self._dir = dirobj

    def set_installed(self):
        self._installed = True

    def is_installed(self):
        return self._installed

    def is_source_file(self):
        return self._sourcefile

    def is_test_file(self):
        return self._dir.is_test_directory()

    def is_documented(self):
        return self._rawdoc.is_documented()

    def get_path(self):
        return self._rawdoc.get_path()

    def get_documentation_type(self):
        return self._rawdoc.get_visibility()

    def get_api_type(self):
        return self._apitype

class Directory(object):
    def __init__(self, rawdoc):
        self._rawdoc = rawdoc
        self._parent = None
        self._is_test_dir = False
        path = rawdoc.get_path().rstrip('/')
        if os.path.basename(path) in ('tests', 'legacytests'):
            self._is_test_dir = True
        self._subdirs = set()

    def _set_parent(self, parent):
        self._parent = parent
        self._is_test_dir = parent.is_test_directory()

    def is_test_directory(self):
        return self._is_test_dir

    def process_subdirectories(self, dirdict):
        for dirdoc in self._rawdoc.get_subdirectories():
            path = dirdoc.get_path()
            dirobj = Directory(dirdoc)
            dirobj._set_parent(self)
            dirdict[path] = dirobj
            self._subdirs.add(dirobj)
            dirobj.process_subdirectories(dirdict)

class SourceTree(object):
    def __init__(self, source_root, build_root, reporter):
        xmldir = os.path.join(build_root, 'doxygen', 'xml')
        docset = xml.DocumentationSet(xmldir)
        docset.load_details()
        docset.remove_duplicates()
        self._docset = docset
        self._dirs = dict()
        rootdirs = docset.get_compounds(xml.Directory, lambda x: x.get_parent() is None)
        for dirdoc in rootdirs:
            path = dirdoc.get_path()
            dirobj = Directory(dirdoc)
            self._dirs[path] = dirobj
            dirobj.process_subdirectories(self._dirs)
        self._files = dict()
        for filedoc in docset.get_files():
            path = filedoc.get_path()
            if not os.path.isabs(path):
                reporter.xml_assert(filedoc.get_xml_path(), "expected absolute path in Doxygen-produced XML vile")
            fileobj = File(filedoc)
            fileobj.set_directory(self._dirs[os.path.dirname(path) + '/'])
            self._files[path] = fileobj

    def set_installed_header_list(self, installedheaders, reporter):
        for path in installedheaders:
            if not os.path.isabs(path):
                reporter.input_error("installed headers should be specified with absolute path: {0}".format(path))
                continue
            if path not in self._files:
                # TODO: After the script handles also files excluded from Doxygen
                # (for include dependency checks), this check makes sense.
                #reporter.input_error("installed header not known to Doxygen: {0}".format(path))
                continue
            self._files[path].set_installed()

    def get_files(self):
        return self._files.itervalues()
