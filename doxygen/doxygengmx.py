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
    def __init__(self, rawdoc, directory, sourcetree):
        self._rawdoc = rawdoc
        self._installed = False
        extension = os.path.splitext(rawdoc.get_path())[1]
        self._sourcefile = False
        if extension in ('.c', '.cpp', '.cu'):
            self._sourcefile = True
        self._dir = directory
        self._apitype = DocType.none
        self._modules = set()
        if self._rawdoc.is_documented():
            grouplist = self._rawdoc.get_groups()
            self._apitype = _get_api_type_for_compound(grouplist)
            for group in grouplist:
                module = sourcetree.get_object(group)
                if module:
                    self._modules.add(module)

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

    def has_brief_description(self):
        return self._rawdoc.has_brief_description()

    def get_path(self):
        return self._rawdoc.get_path()

    def get_documentation_type(self):
        return self._rawdoc.get_visibility()

    def get_api_type(self):
        return self._apitype

    def get_expected_module(self):
        return self._dir.get_module()

    def get_doc_modules(self):
        return self._modules

class Directory(object):
    def __init__(self, rawdoc):
        self._rawdoc = rawdoc
        self._parent = None
        self._module = None
        self._is_test_dir = False
        path = rawdoc.get_path().rstrip('/')
        if os.path.basename(path) in ('tests', 'legacytests'):
            self._is_test_dir = True
        self._subdirs = set()

    def _set_parent(self, parent):
        self._parent = parent
        self._is_test_dir = parent.is_test_directory()

    def set_module(self, module):
        self._module = module

    def is_test_directory(self):
        return self._is_test_dir

    def get_module(self):
        if self._module:
            return self._module
        if self._parent:
            return self._parent.get_module()
        return None

    def process_subdirectories(self, dirdict, objmap):
        for dirdoc in self._rawdoc.get_subdirectories():
            path = dirdoc.get_path()
            dirobj = Directory(dirdoc)
            dirobj._set_parent(self)
            objmap[dirdoc] = dirobj
            dirdict[path] = dirobj
            self._subdirs.add(dirobj)
            dirobj.process_subdirectories(dirdict, objmap)

class Module(object):
    def __init__(self, rawdoc, rootdir):
        self._rawdoc = rawdoc
        self._rootdir = rootdir

    def get_name(self):
        return self._rawdoc.get_name()

class SourceTree(object):
    def __init__(self, source_root, build_root, reporter):
        self._source_root = source_root
        self._reporter = reporter
        xmldir = os.path.join(build_root, 'doxygen', 'xml')
        docset = xml.DocumentationSet(xmldir)
        docset.load_details()
        docset.remove_duplicates()
        self._docset = docset
        self._objmap = dict()
        self._load_dirs()
        self._load_modules()
        self._files = dict()
        for filedoc in docset.get_files():
            path = filedoc.get_path()
            if not os.path.isabs(path):
                reporter.xml_assert(filedoc.get_xml_path(), "expected absolute path in Doxygen-produced XML file")
            extension = os.path.splitext(filedoc.get_path())[1]
            if extension == '.md':
                continue
            dirdoc = filedoc.get_directory()
            if not dirdoc:
                reporter.xml_assert(filedoc.get_xml_path(), "file is not in any directory in Doxygen")
                continue
            fileobj = File(filedoc, self._objmap[dirdoc], self)
            self._objmap[filedoc] = fileobj
            self._files[path] = fileobj

    def _load_dirs(self):
        self._dirs = dict()
        rootdirs = self._docset.get_compounds(xml.Directory, lambda x: x.get_parent() is None)
        for dirdoc in rootdirs:
            path = dirdoc.get_path()
            dirobj = Directory(dirdoc)
            self._objmap[dirdoc] = dirobj
            self._dirs[path] = dirobj
            dirobj.process_subdirectories(self._dirs, self._objmap)

    def _load_modules(self):
        self._modules = set()
        moduledocs = self._docset.get_compounds(xml.Group,
                lambda x: x.get_name().startswith('module_'))
        for moduledoc in moduledocs:
            # Remove module_ prefix
            basename = moduledoc.get_name()[7:]
            if '_' in basename:
                # TODO: Consider a better rule
                continue
            rootdir = self._get_dir(os.path.join('src', 'gromacs', basename))
            if not rootdir:
                rootdir = self._get_dir(os.path.join('src', basename))
            if not rootdir:
                self._reporter.input_error("no matching directory for module: {0}".format(moduledoc))
                continue
            moduleobj = Module(moduledoc, rootdir)
            rootdir.set_module(moduleobj)
            self._objmap[moduledoc] = moduleobj

    def _get_dir(self, relpath):
        return self._dirs.get(os.path.join(self._source_root, relpath, ''))

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

    def get_object(self, docobj):
        return self._objmap.get(docobj)

    def get_files(self):
        return self._files.itervalues()

    def get_code_entities(self):
        result = set(self._docset.get_members())
        result.update(self._docset.get_compounds((xml.Namespace, xml.Class)))
        return result
