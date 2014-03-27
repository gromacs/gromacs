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

"""GROMACS-specific representation for source tree and documentation.

This module provides classes that construct a GROMACS-specific representation
of the source tree and associate the Doxygen XML output with it.  It constructs
an initial representation by walking the source tree in the file system, and
then associates information from the Doxygen XML output into this.
It also adds some additional knowledge from how the GROMACS source tree is
organized to construct a representation that is easy to process and check as
the top-level scripts expect.

The object model is rooted at a GromacsTree object.  Currently, it constructs a
representation of the source tree from the file system, but is otherwise mostly
a thin wrapper around the Doxygen XML tree.  It already adds some relations and
rules that come from GROMACS-specific knowledge.  In the future, more such
customizations will be added.
"""

import os
import os.path

import doxygenxml as xml
import reporter
# We import DocType directly so that it is exposed from this module as well.
from doxygenxml import DocType

def _get_api_type_for_compound(grouplist):
    """Helper function to deduce API type from Doxygen group membership."""
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

    """Source/header file in the GROMACS tree."""

    def __init__(self, path, directory):
        """Initialize a file representation with basic information."""
        self._path = path
        self._dir = directory
        self._rawdoc = None
        self._installed = False
        extension = os.path.splitext(path)[1]
        self._sourcefile = (extension in ('.c', '.cc', '.cpp', '.cu'))
        self._apitype = DocType.none
        self._modules = set()

    def set_doc_xml(self, rawdoc, sourcetree):
        """Assiociate Doxygen documentation entity with the file."""
        assert self._rawdoc is None
        assert rawdoc.is_source_file() == self._sourcefile
        self._rawdoc = rawdoc
        if self._rawdoc.is_documented():
            grouplist = self._rawdoc.get_groups()
            self._apitype = _get_api_type_for_compound(grouplist)
            for group in grouplist:
                module = sourcetree.get_object(group)
                if module:
                    self._modules.add(module)

    def set_installed(self):
        """Mark the file installed."""
        self._installed = True

    def get_reporter_location(self):
        return reporter.Location(self._path, None)

    def is_installed(self):
        return self._installed

    def is_source_file(self):
        return self._sourcefile

    def is_test_file(self):
        return self._dir.is_test_directory()

    def is_documented(self):
        return self._rawdoc and self._rawdoc.is_documented()

    def has_brief_description(self):
        return self._rawdoc and self._rawdoc.has_brief_description()

    def get_path(self):
        return self._path

    def get_documentation_type(self):
        if not self._rawdoc:
            return DocType.none
        return self._rawdoc.get_visibility()

    def get_api_type(self):
        return self._apitype

    def get_expected_module(self):
        return self._dir.get_module()

    def get_doc_modules(self):
        return self._modules

class GeneratedFile(File):
    pass

class Directory(object):

    """(Sub)directory in the GROMACS tree."""

    def __init__(self, path, parent):
        """Initialize a file representation with basic information."""
        self._path = path
        self._name = os.path.basename(path)
        self._parent = parent
        self._rawdoc = None
        self._module = None
        self._is_test_dir = False
        if parent and parent.is_test_directory() or \
                os.path.basename(path) in ('tests', 'legacytests'):
            self._is_test_dir = True
        self._subdirs = set()
        if parent:
            parent._subdirs.add(self)

    def set_doc_xml(self, rawdoc, sourcetree):
        """Assiociate Doxygen documentation entity with the directory."""
        assert self._rawdoc is None
        assert self._path == rawdoc.get_path().rstrip('/')
        self._rawdoc = rawdoc

    def set_module(self, module):
        assert self._module is None
        self._module = module

    def get_name(self):
        return self._name

    def get_reporter_location(self):
        return reporter.Location(self._path, None)

    def is_test_directory(self):
        return self._is_test_dir

    def get_module(self):
        if self._module:
            return self._module
        if self._parent:
            return self._parent.get_module()
        return None

    def get_subdirectories(self):
        return self._subdirs

class Module(object):

    """Code module in the GROMACS source tree.

    Modules are specific subdirectories that host a more or less coherent
    set of routines.  Simplified, every subdirectory under src/gromacs/ is
    a different module.  This object provides that abstraction and also links
    the subdirectory to the module documentation (documented as a group in
    Doxygen) if that exists.
    """

    def __init__(self, name, rootdir):
        self._name = name
        self._rawdoc = None
        self._rootdir = rootdir

    def set_doc_xml(self, rawdoc, sourcetree):
        """Assiociate Doxygen documentation entity with the module."""
        assert self._rawdoc is None
        self._rawdoc = rawdoc

    def is_documented(self):
        return self._rawdoc is not None

    def get_name(self):
        return self._name

class Class(object):

    """Class/struct/union in the GROMACS source code."""

    def __init__(self, rawdoc, files):
        self._rawdoc = rawdoc
        self._files = set(files)

    def get_name(self):
        return self._rawdoc.get_name()

    def get_reporter_location(self):
        return self._rawdoc.get_reporter_location()

    def get_files(self):
        return self._files

    def is_documented(self):
        return self._rawdoc.is_documented()

    def has_brief_description(self):
        return self._rawdoc.has_brief_description()

    def get_documentation_type(self):
        if not self.is_documented():
            return DocType.none
        if self._rawdoc.is_local():
            return DocType.internal
        return self._rawdoc.get_visibility()

    def get_file_documentation_type(self):
        return max([fileobj.get_documentation_type() for fileobj in self._files])

    def is_in_installed_file(self):
        return any([fileobj.is_installed() for fileobj in self._files])

class GromacsTree(object):

    """Root object for navigating the GROMACS source tree.

    On initialization, the list of files and directories is initialized by
    walking the source tree, and modules are created for top-level
    subdirectories.  At this point, only information that is accessible from
    file names and paths only is available.

    set_installed_file_list() can be called to set the list of installed
    files.

    load_xml() can be called to load information from Doxygen XML data in
    the build tree (the Doxygen XML data must have been built separately).
    """

    def __init__(self, source_root, build_root, reporter):
        """Initialize the tree object by walking the source tree."""
        self._source_root = os.path.abspath(source_root)
        self._build_root = os.path.abspath(build_root)
        self._reporter = reporter
        self._docset = None
        self._docmap = dict()
        self._dirs = dict()
        self._files = dict()
        self._modules = dict()
        self._classes = set()
        self._walk_dir(os.path.join(self._source_root, 'src'))
        rootdir = self._get_dir(os.path.join('src', 'gromacs'))
        for subdir in rootdir.get_subdirectories():
            self._create_module(subdir)
        rootdir = self._get_dir(os.path.join('src', 'testutils'))
        self._create_module(rootdir)

    def _get_rel_path(self, path):
        assert os.path.isabs(path)
        if path.startswith(self._build_root):
            return path[len(self._build_root)+1:]
        if path.startswith(self._source_root):
            return path[len(self._source_root)+1:]
        raise ValueError("path not under build nor source tree: {0}".format(path))

    def _walk_dir(self, rootpath):
        """Construct representation of the source tree by walking the file system."""
        assert os.path.isabs(rootpath)
        assert rootpath not in self._dirs
        relpath = self._get_rel_path(rootpath)
        self._dirs[relpath] = Directory(rootpath, None)
        for dirpath, dirnames, filenames in os.walk(rootpath):
            if 'contrib' in dirnames:
                dirnames.remove('contrib')
            if 'refdata' in dirnames:
                dirnames.remove('refdata')
            currentdir = self._dirs[self._get_rel_path(dirpath)]
            # Loop through a copy so that we can modify dirnames.
            for dirname in list(dirnames):
                fullpath = os.path.join(dirpath, dirname)
                if fullpath == self._build_root:
                    dirnames.remove(dirname)
                    continue
                relpath = self._get_rel_path(fullpath)
                self._dirs[relpath] = Directory(fullpath, currentdir)
            extensions = ('.h', '.cuh', '.hpp', '.c', '.cc', '.cpp', '.cu')
            for filename in filenames:
                basename, extension = os.path.splitext(filename)
                if extension in extensions:
                    fullpath = os.path.join(dirpath, filename)
                    relpath = self._get_rel_path(fullpath)
                    self._files[relpath] = File(fullpath, currentdir)
                elif extension == '.cmakein':
                    extension = os.path.splitext(basename)[1]
                    if extension in extensions:
                        fullpath = os.path.join(dirpath, basename)
                        relpath = self._get_rel_path(fullpath)
                        fullpath = os.path.join(dirpath, filename)
                        self._files[relpath] = GeneratedFile(fullpath, currentdir)

    def _create_module(self, rootdir):
        """Create module for a subdirectory."""
        name = 'module_' + rootdir.get_name()
        moduleobj = Module(name, rootdir)
        rootdir.set_module(moduleobj)
        self._modules[name] = moduleobj

    def load_xml(self):
        """Load Doxygen XML information."""
        xmldir = os.path.join(self._build_root, 'doxygen', 'xml')
        self._docset = xml.DocumentationSet(xmldir, self._reporter)
        self._docset.load_details()
        self._docset.merge_duplicates()
        self._load_dirs()
        self._load_modules()
        self._load_files()
        self._load_classes()

    def _load_dirs(self):
        """Load Doxygen XML directory information."""
        rootdirs = self._docset.get_compounds(xml.Directory,
                lambda x: x.get_parent() is None)
        for dirdoc in rootdirs:
            self._load_dir(dirdoc, None)

    def _load_dir(self, dirdoc, parent):
        """Load Doxygen XML directory information for a single directory."""
        path = dirdoc.get_path().rstrip('/')
        if not os.path.isabs(path):
            self._reporter.xml_assert(dirdoc.get_xml_path(),
                    "expected absolute path in Doxygen-produced XML file")
            return
        relpath = self._get_rel_path(path)
        dirobj = self._dirs.get(relpath)
        if not dirobj:
            dirobj = Directory(path, parent)
            self._dirs[relpath] = dirobj
        dirobj.set_doc_xml(dirdoc, self)
        self._docmap[dirdoc] = dirobj
        for subdirdoc in dirdoc.get_subdirectories():
            self._load_dir(subdirdoc, dirobj)

    def _load_modules(self):
        """Load Doxygen XML module (group) information."""
        moduledocs = self._docset.get_compounds(xml.Group,
                lambda x: x.get_name().startswith('module_'))
        for moduledoc in moduledocs:
            moduleobj = self._modules.get(moduledoc.get_name())
            if not moduleobj:
                self._reporter.input_error(
                        "no matching directory for module: {0}".format(moduledoc))
                continue
            moduleobj.set_doc_xml(moduledoc, self)
            self._docmap[moduledoc] = moduleobj

    def _load_files(self):
        """Load Doxygen XML file information."""
        for filedoc in self._docset.get_files():
            path = filedoc.get_path()
            if not os.path.isabs(path):
                self._reporter.xml_assert(filedoc.get_xml_path(),
                        "expected absolute path in Doxygen-produced XML file")
                continue
            extension = os.path.splitext(filedoc.get_path())[1]
            # We don't care about Markdown files that only produce pages
            # (and fail the directory check below).
            if extension == '.md':
                continue
            dirdoc = filedoc.get_directory()
            if not dirdoc:
                self._reporter.xml_assert(filedoc.get_xml_path(),
                        "file is not in any directory in Doxygen")
                continue
            relpath = self._get_rel_path(path)
            fileobj = self._files.get(relpath)
            if not fileobj:
                fileobj = File(path, self._docmap[dirdoc])
                self._files[relpath] = fileobj
            fileobj.set_doc_xml(filedoc, self)
            self._docmap[filedoc] = fileobj

    def _load_classes(self):
        """Load Doxygen XML class information."""
        classdocs = self._docset.get_classes()
        for classdoc in classdocs:
            files = [self._docmap[filedoc] for filedoc in classdoc.get_files()]
            classobj = Class(classdoc, files)
            self._docmap[classdoc] = classobj
            self._classes.add(classobj)

    def _get_dir(self, relpath):
        """Get directory object for a path relative to source tree root."""
        return self._dirs.get(relpath)

    def set_installed_file_list(self, installedfiles):
        """Set list of installed files."""
        for path in installedfiles:
            if not os.path.isabs(path):
                self._reporter.input_error(
                        "installed file not specified with absolute path: {0}"
                        .format(path))
                continue
            relpath = self._get_rel_path(path)
            if relpath not in self._files:
                self._reporter.input_error(
                        "installed file not in source tree: {0}".format(path))
                continue
            self._files[relpath].set_installed()

    def get_object(self, docobj):
        """Get tree object for a Doxygen XML object."""
        return self._docmap.get(docobj)

    def get_files(self):
        """Get iterable for all files in the source tree."""
        return self._files.itervalues()

    def get_classes(self):
        """Get iterable for all classes in the source tree."""
        return self._classes

    def get_members(self):
        """Get iterable for all members (in Doxygen terms) in the source tree."""
        # TODO: Add wrappers to solve some issues.
        return self._docset.get_members()
