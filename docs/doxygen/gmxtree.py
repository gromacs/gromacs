#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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

import collections
import os
import os.path
import re
import subprocess

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

class IncludedFile(object):

    """Information about an #include directive in a file."""

    def __init__(self, including_file, lineno, included_file, included_path, is_relative, is_system, line):
        self._including_file = including_file
        self._line_number = lineno
        self._included_file = included_file
        self._included_path = included_path
        #self._used_include_path = used_include_path
        self._is_relative = is_relative
        self._is_system = is_system
        self._line = line

    def __str__(self):
        if self._is_system:
            return '<{0}>'.format(self._included_path)
        else:
            return '"{0}"'.format(self._included_path)

    def is_system(self):
        return self._is_system

    def is_relative(self):
        return self._is_relative

    def get_included_path(self):
        return self._included_path

    def get_including_file(self):
        return self._including_file

    def get_file(self):
        return self._included_file

    def get_line_number(self):
        return self._line_number

    def get_full_line(self):
        """Return the full source line on which this include appears.

        Trailing newline is included."""
        return self._line

    def get_reporter_location(self):
        return reporter.Location(self._including_file.get_abspath(), self._line_number)

class IncludeBlock(object):

    """Block of consequent #include directives in a file."""

    def __init__(self, first_included_file):
        self._first_line = first_included_file.get_line_number()
        self._last_line = self._first_line
        self._files = []
        self.add_file(first_included_file)

    def add_file(self, included_file):
        self._files.append(included_file)
        self._last_line = included_file.get_line_number()

    def get_includes(self):
        return self._files

    def get_first_line(self):
        return self._first_line

    def get_last_line(self):
        return self._last_line

class File(object):

    """Source/header file in the GROMACS tree."""

    def __init__(self, abspath, relpath, directory):
        """Initialize a file representation with basic information."""
        self._abspath = abspath
        self._relpath = relpath
        self._dir = directory
        self._rawdoc = None
        self._installed = False
        extension = os.path.splitext(abspath)[1]
        self._sourcefile = (extension in ('.c', '.cc', '.cpp', '.cu'))
        self._apitype = DocType.none
        self._modules = set()
        self._includes = []
        self._include_blocks = []
        self._main_header = None
        self._lines = None
        self._filter = None
        self._declared_defines = None
        self._used_define_files = set()
        directory.add_file(self)

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

    def set_git_filter_attribute(self, filtername):
        """Set the git filter attribute associated with the file."""
        self._filter = filtername

    def set_main_header(self, included_file):
        """Set the main header file for a source file."""
        assert self.is_source_file()
        self._main_header = included_file

    def _process_include(self, lineno, is_system, includedpath, line, sourcetree):
        """Process #include directive during scan()."""
        is_relative = False
        if is_system:
            fileobj = sourcetree.find_include_file(includedpath)
        else:
            fullpath = os.path.join(self._dir.get_abspath(), includedpath)
            fullpath = os.path.abspath(fullpath)
            if os.path.exists(fullpath):
                is_relative = True
                fileobj = sourcetree.get_file(fullpath)
            else:
                fileobj = sourcetree.find_include_file(includedpath)
        included_file = IncludedFile(self, lineno, fileobj, includedpath,
            is_relative, is_system, line)
        self._includes.append(included_file)
        return included_file

    def scan_contents(self, sourcetree, keep_contents, detect_defines):
        """Scan the file contents and initialize information based on it."""
        # TODO: Consider a more robust regex.
        include_re = r'^\s*#\s*include\s+(?P<quote>["<])(?P<path>[^">]*)[">]'
        define_re = r'^\s*#.*define\s+(\w*)'
        current_block = None
        with open(self._abspath, 'r') as scanfile:
            contents = scanfile.read()
        lines = contents.splitlines(True)
        for lineno, line in enumerate(lines, 1):
            match = re.match(include_re, line)
            if match:
                is_system = (match.group('quote') == '<')
                includedpath = match.group('path')
                included_file = self._process_include(lineno, is_system,
                        includedpath, line, sourcetree)
                if current_block is None:
                    current_block = IncludeBlock(included_file)
                    self._include_blocks.append(current_block)
                else:
                    current_block.add_file(included_file)
            elif line and not line.isspace():
                current_block = None
        if detect_defines:
            self._declared_defines = []
            for line in lines:
                match = re.match(define_re, line)
                if match:
                    self._declared_defines.append(match.group(1))
        if keep_contents:
            self._lines = lines

    def add_used_defines(self, define_file, defines):
        """Add defines used in this file.

        Used internally by find_define_file_uses()."""
        self._used_define_files.add(define_file)

    def get_reporter_location(self):
        return reporter.Location(self._abspath, None)

    def is_installed(self):
        return self._installed

    def is_external(self):
        return self._dir.is_external()

    def is_source_file(self):
        return self._sourcefile

    def is_test_file(self):
        return self._dir.is_test_directory()

    def should_includes_be_sorted(self):
        """Return whether the include directives in the file should be sorted."""
        return self._filter in ('includesort', 'uncrustify')

    def is_documented(self):
        return self._rawdoc and self._rawdoc.is_documented()

    def has_brief_description(self):
        return self._rawdoc and self._rawdoc.has_brief_description()

    def get_abspath(self):
        return self._abspath

    def get_relpath(self):
        return self._relpath

    def get_name(self):
        return os.path.basename(self._abspath)

    def get_directory(self):
        return self._dir

    def get_doc_type(self):
        if not self._rawdoc:
            return DocType.none
        return self._rawdoc.get_visibility()

    def get_api_type(self):
        return self._apitype

    def api_type_is_reliable(self):
        if self._apitype in (DocType.internal, DocType.library):
            return True
        module = self.get_module()
        return module and module.is_documented()

    def is_public(self):
        if self.api_type_is_reliable():
            return self.get_api_type() == DocType.public
        return self.get_api_type() == DocType.public or self.is_installed()

    def is_module_internal(self):
        if self.is_source_file():
            return True
        return not self.is_installed() and self.get_api_type() <= DocType.internal

    def get_expected_module(self):
        return self._dir.get_module()

    def get_doc_modules(self):
        return self._modules

    def get_module(self):
        module = self.get_expected_module()
        if not module and len(self._modules) == 1:
            module = list(self._modules)[0]
        return module

    def get_includes(self):
        return self._includes

    def get_include_blocks(self):
        return self._include_blocks

    def _get_included_files_recurse(self, result):
        for include in self._includes:
            included_file = include.get_file()
            if included_file is not None and not included_file in result:
                result.add(included_file)
                included_file._get_included_files_recurse(result)

    def get_included_files(self, recursive=False):
        if recursive:
            result = set()
            self._get_included_files_recurse(result)
            return result
        return set([x.get_file() for x in self._includes])

    def get_main_header(self):
        return self._main_header

    def get_contents(self):
        return self._lines

    def get_declared_defines(self):
        """Return set of defines declared in this file.

        The information is only populated for selected files."""
        return self._declared_defines

    def get_used_define_files(self):
        """Return set of defines from config.h that are used in this file.

        The return value is empty if find_define_file_uses() has not been called,
        as well as for headers that declare these defines."""
        return self._used_define_files

class GeneratedFile(File):
    def __init__(self, abspath, relpath, directory):
        File.__init__(self, abspath, relpath, directory)
        self._generator_source_file = None

    def scan_contents(self, sourcetree, keep_contents, detect_defines):
        if os.path.exists(self.get_abspath()):
            File.scan_contents(self, sourcetree, keep_contents, False)

    def set_generator_source(self, sourcefile):
        self._generator_source_file = sourcefile

    def get_generator_source(self):
        return self._generator_source_file

    def get_reporter_location(self):
        if self._generator_source_file:
            return self._generator_source_file.get_reporter_location()
        return File.get_reporter_location(self)

    def get_declared_defines(self):
        if self._generator_source_file:
            return self._generator_source_file.get_declared_defines()
        return File.get_declared_defines(self)

class GeneratorSourceFile(File):
    pass

class Directory(object):

    """(Sub)directory in the GROMACS tree."""

    def __init__(self, abspath, relpath, parent):
        """Initialize a file representation with basic information."""
        self._abspath = abspath
        self._relpath = relpath
        self._name = os.path.basename(abspath)
        self._parent = parent
        self._rawdoc = None
        self._module = None
        self._is_test_dir = False
        if parent and parent.is_test_directory() or \
                self._name in ('tests', 'legacytests'):
            self._is_test_dir = True
        self._is_external = False
        if parent and parent.is_external() or self._name == 'external':
            self._is_external = True
        self._subdirs = set()
        if parent:
            parent._subdirs.add(self)
        self._files = set()
        self._has_installed_files = None

    def set_doc_xml(self, rawdoc, sourcetree):
        """Assiociate Doxygen documentation entity with the directory."""
        assert self._rawdoc is None
        assert self._abspath == rawdoc.get_path().rstrip('/')
        self._rawdoc = rawdoc

    def set_module(self, module):
        assert self._module is None
        self._module = module

    def add_file(self, fileobj):
        self._files.add(fileobj)

    def get_name(self):
        return self._name

    def get_reporter_location(self):
        return reporter.Location(self._abspath, None)

    def get_abspath(self):
        return self._abspath

    def get_relpath(self):
        return self._relpath

    def is_test_directory(self):
        return self._is_test_dir

    def is_external(self):
        return self._is_external

    def has_installed_files(self):
        if self._has_installed_files is None:
            self._has_installed_files = False
            for subdir in self._subdirs:
                if subdir.has_installed_files():
                    self._has_installed_files = True
                    return True
            for fileobj in self._files:
                if fileobj.is_installed():
                    self._has_installed_files = True
                    return True
        return self._has_installed_files

    def get_module(self):
        if self._module:
            return self._module
        if self._parent:
            return self._parent.get_module()
        return None

    def get_subdirectories(self):
        return self._subdirs

    def get_files(self):
        for subdir in self._subdirs:
            for fileobj in subdir.get_files():
                yield fileobj
        for fileobj in self._files:
            yield fileobj

    def contains(self, fileobj):
        """Check whether file is within the directory or its subdirectories."""
        dirobj = fileobj.get_directory()
        while dirobj:
            if dirobj == self:
                return True
            dirobj = dirobj._parent
        return False

class ModuleDependency(object):

    """Dependency between modules."""

    def __init__(self, othermodule):
        """Initialize empty dependency object with given module as dependency."""
        self._othermodule = othermodule
        self._includedfiles = []
        self._cyclesuppression = None
        self._is_test_only_dependency = True

    def add_included_file(self, includedfile):
        """Add IncludedFile that is part of this dependency."""
        assert includedfile.get_file().get_module() == self._othermodule
        if not includedfile.get_including_file().is_test_file():
            self._is_test_only_dependency = False
        self._includedfiles.append(includedfile)

    def set_cycle_suppression(self):
        """Set suppression on cycles containing this dependency."""
        self._cyclesuppression = True

    def is_cycle_suppressed(self):
        """Return whether cycles containing this dependency are suppressed."""
        return self._cyclesuppression is not None

    def is_test_only_dependency(self):
        """Return whether this dependency is only from test code."""
        return self._is_test_only_dependency

    def get_other_module(self):
        """Get module that this dependency is to."""
        return self._othermodule

    def get_included_files(self):
        """Get IncludedFile objects for the individual include dependencies."""
        return self._includedfiles

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
        self._group = None
        self._dependencies = dict()

    def set_doc_xml(self, rawdoc, sourcetree):
        """Assiociate Doxygen documentation entity with the module."""
        assert self._rawdoc is None
        self._rawdoc = rawdoc
        if self._rawdoc.is_documented():
            groups = list(self._rawdoc.get_groups())
            if len(groups) == 1:
                groupname = groups[0].get_name()
                if groupname.startswith('group_'):
                    self._group = groupname[6:]

    def add_dependency(self, othermodule, includedfile):
        """Add #include dependency from a file in this module."""
        assert includedfile.get_file().get_module() == othermodule
        if othermodule not in self._dependencies:
            self._dependencies[othermodule] = ModuleDependency(othermodule)
        self._dependencies[othermodule].add_included_file(includedfile)

    def is_documented(self):
        return self._rawdoc is not None

    def get_name(self):
        return self._name

    def get_root_dir(self):
        return self._rootdir

    def get_files(self):
        # TODO: Include public API convenience headers?
        return self._rootdir.get_files()

    def get_group(self):
        return self._group

    def get_dependencies(self):
        return self._dependencies.itervalues()

class Namespace(object):

    """Namespace in the GROMACS source code."""

    def __init__(self, rawdoc):
        self._rawdoc = rawdoc

    def is_anonymous(self):
        return self._rawdoc.is_anonymous()

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

    def get_doc_type(self):
        """Return documentation type (visibility) for the class.

        In addition to the actual code, this encodes GROMACS-specific logic
        of setting EXTRACT_LOCAL_CLASSES=YES only for the full documentation.
        Local classes never appear outside the full documentation, no matter
        what is their visibility.
        """
        if not self.is_documented():
            return DocType.none
        if self._rawdoc.is_local():
            return DocType.internal
        return self._rawdoc.get_visibility()

    def get_file_doc_type(self):
        return max([fileobj.get_doc_type() for fileobj in self._files])

    def is_in_installed_file(self):
        return any([fileobj.is_installed() for fileobj in self._files])

class Member(object):

    """Member (in Doxygen terminology) in the GROMACS source tree.

    Currently, modeling is limited to the minimal set of properties that the
    checker uses.
    """

    def __init__(self, rawdoc, namespace):
        self._rawdoc = rawdoc
        self._namespace = namespace

    def get_name(self):
        return self._rawdoc.get_name()

    def get_reporter_location(self):
        return self._rawdoc.get_reporter_location()

    def is_documented(self):
        return self._rawdoc.is_documented()

    def has_brief_description(self):
        return self._rawdoc.has_brief_description()

    def has_inbody_description(self):
        return self._rawdoc.has_inbody_description()

    def is_visible(self):
        """Return whether the member is visible in Doxygen documentation.

        Doxygen ignores members whose parent compounds are not documented.
        However, when EXTRACT_ANON_NPACES=ON (which is set for our full
        documentation), members of anonymous namespaces are extracted even if
        the namespace is the only parent and is not documented.
        """
        if self._namespace and self._namespace.is_anonymous():
            return True
        return self._rawdoc.get_inherited_visibility() != DocType.none


class GromacsTree(object):

    """Root object for navigating the GROMACS source tree.

    On initialization, the list of files and directories is initialized by
    walking the source tree, and modules are created for top-level
    subdirectories.  At this point, only information that is accessible from
    file names and paths only is available.

    load_git_attributes() can be called to load attribute information from
    .gitattributes for all the files.

    load_installed_file_list() can be called to load the list of installed
    files from the build tree (generated by CMake).

    scan_files() can be called to read all the files and initialize #include
    dependencies between the files based on the information.  This is done like
    this instead of relying on Doxygen-extracted include files to make the
    dependency graph independent from preprocessor macro definitions
    (Doxygen only sees those #includes that the preprocessor sees, which
    depends on what #defines it has seen).

    find_define_file_uses() can be called to find all uses of defines
    declared in config.h and some other macro headers. In the current
    implementation, scan_files() must have been called earlier.

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
        self._namespaces = set()
        self._members = set()
        self._walk_dir(os.path.join(self._source_root, 'src'))
        for fileobj in self.get_files():
            if fileobj and fileobj.is_source_file() and not fileobj.is_external():
                (basedir, name) = os.path.split(fileobj.get_abspath())
                (basename, ext) = os.path.splitext(name)
                header = self.get_file(os.path.join(basedir, basename + '.h'))
                if not header and ext == '.cu':
                    header = self.get_file(os.path.join(basedir, basename + '.cuh'))
                if not header and fileobj.is_test_file():
                    basedir = os.path.dirname(basedir)
                    header = self.get_file(os.path.join(basedir, basename + '.h'))
                    if not header:
                        # Somewhat of a hack; currently, the tests for
                        # analysisdata/modules/ and trajectoryanalysis/modules/
                        # is at the top-level tests directory.
                        # TODO: It could be clearer to split the tests so that
                        # there would be a separate modules/tests/.
                        header = self.get_file(os.path.join(basedir, 'modules', basename + '.h'))
                    if not header and basename.endswith('_tests'):
                        header = self.get_file(os.path.join(basedir, basename[:-6] + '.h'))
                if not header and fileobj.get_relpath().startswith('src/gromacs'):
                    header = self._files.get(os.path.join('src/gromacs/legacyheaders', basename + '.h'))
                if header:
                    fileobj.set_main_header(header)
        rootdir = self._get_dir(os.path.join('src', 'gromacs'))
        for subdir in rootdir.get_subdirectories():
            self._create_module(subdir)
        rootdir = self._get_dir(os.path.join('src', 'testutils'))
        self._create_module(rootdir)

    def _get_rel_path(self, path):
        assert os.path.isabs(path)
        if path.startswith(self._build_root):
            return os.path.relpath(path, self._build_root)
        if path.startswith(self._source_root):
            return os.path.relpath(path, self._source_root)
        raise ValueError("path not under build nor source tree: {0}".format(path))

    def _walk_dir(self, rootpath):
        """Construct representation of the source tree by walking the file system."""
        assert os.path.isabs(rootpath)
        assert rootpath not in self._dirs
        relpath = self._get_rel_path(rootpath)
        self._dirs[relpath] = Directory(rootpath, relpath, None)
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
                self._dirs[relpath] = Directory(fullpath, relpath, currentdir)
            extensions = ('.h', '.cuh', '.hpp', '.c', '.cc', '.cpp', '.cu', '.bm')
            for filename in filenames:
                basename, extension = os.path.splitext(filename)
                if extension in extensions:
                    fullpath = os.path.join(dirpath, filename)
                    relpath = self._get_rel_path(fullpath)
                    self._files[relpath] = File(fullpath, relpath, currentdir)
                elif extension == '.cmakein':
                    extension = os.path.splitext(basename)[1]
                    if extension in extensions:
                        fullpath = os.path.join(dirpath, filename)
                        relpath = self._get_rel_path(fullpath)
                        sourcefile = GeneratorSourceFile(fullpath, relpath, currentdir)
                        self._files[relpath] = sourcefile
                        fullpath = os.path.join(dirpath, basename)
                        relpath = self._get_rel_path(fullpath)
                        fullpath = os.path.join(self._build_root, relpath)
                        generatedfile = GeneratedFile(fullpath, relpath, currentdir)
                        self._files[relpath] = generatedfile
                        generatedfile.set_generator_source(sourcefile)
                elif extension in ('.l', '.y', '.pre'):
                    fullpath = os.path.join(dirpath, filename)
                    relpath = self._get_rel_path(fullpath)
                    self._files[relpath] = GeneratorSourceFile(fullpath, relpath, currentdir)

    def _create_module(self, rootdir):
        """Create module for a subdirectory."""
        name = 'module_' + rootdir.get_name()
        moduleobj = Module(name, rootdir)
        rootdir.set_module(moduleobj)
        self._modules[name] = moduleobj

    def scan_files(self, only_files=None, keep_contents=False):
        """Read source files to initialize #include dependencies."""
        if only_files:
            filelist = only_files
        else:
            filelist = self._files.itervalues()
        define_files = list(self.get_checked_define_files())
        for define_file in list(define_files):
            if isinstance(define_file, GeneratedFile) and \
                    define_file.get_generator_source() is not None:
                define_files.append(define_file.get_generator_source())
        for fileobj in filelist:
            if not fileobj.is_external():
                detect_defines = fileobj in define_files
                fileobj.scan_contents(self, keep_contents, detect_defines)
                module = fileobj.get_module()
                if module:
                    for includedfile in fileobj.get_includes():
                        otherfile = includedfile.get_file()
                        if otherfile:
                            othermodule = otherfile.get_module()
                            if othermodule and othermodule != module:
                                module.add_dependency(othermodule, includedfile)

    def load_xml(self, only_files=None):
        """Load Doxygen XML information.

        If only_files is True, XML data is not loaded for code constructs, but
        only for files, directories, and their potential parents.
        """
        xmldir = os.path.join(self._build_root, 'docs', 'html', 'doxygen', 'xml')
        self._docset = xml.DocumentationSet(xmldir, self._reporter)
        if only_files:
            if isinstance(only_files, collections.Iterable):
                filelist = [x.get_abspath() for x in only_files]
                self._docset.load_file_details(filelist)
            else:
                self._docset.load_file_details()
        else:
            self._docset.load_details()
            self._docset.merge_duplicates()
        self._load_dirs()
        self._load_modules()
        self._load_files()
        if not only_files:
            self._load_namespaces()
            self._load_classes()
            self._load_members()

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
            dirobj = Directory(path, relpath, parent)
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
            if not path:
                # In case of only partially loaded file information,
                # the path information is not set for unloaded files.
                continue
            if not os.path.isabs(path):
                self._reporter.xml_assert(filedoc.get_xml_path(),
                        "expected absolute path in Doxygen-produced XML file")
                continue
            extension = os.path.splitext(path)[1]
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
                fileobj = File(path, relpath, self._docmap[dirdoc])
                self._files[relpath] = fileobj
            fileobj.set_doc_xml(filedoc, self)
            self._docmap[filedoc] = fileobj

    def _load_namespaces(self):
        """Load Doxygen XML namespace information."""
        nsdocs = self._docset.get_namespaces()
        for nsdoc in nsdocs:
            nsobj = Namespace(nsdoc)
            self._docmap[nsdoc] = nsobj
            self._namespaces.add(nsobj)

    def _load_classes(self):
        """Load Doxygen XML class information."""
        classdocs = self._docset.get_classes()
        for classdoc in classdocs:
            files = [self._docmap[filedoc] for filedoc in classdoc.get_files()]
            classobj = Class(classdoc, files)
            self._docmap[classdoc] = classobj
            self._classes.add(classobj)

    def _load_members(self):
        """Load Doxygen XML member information."""
        memberdocs = self._docset.get_members()
        for memberdoc in memberdocs:
            nsdoc = memberdoc.get_namespace()
            nsobj = self.get_object(nsdoc)
            memberobj = Member(memberdoc, nsobj)
            self._docmap[memberdoc] = memberobj
            self._members.add(memberobj)

    def _get_dir(self, relpath):
        """Get directory object for a path relative to source tree root."""
        return self._dirs.get(relpath)

    def get_file(self, path):
        """Get file object for a path relative to source tree root."""
        return self._files.get(self._get_rel_path(path))

    def find_include_file(self, includedpath):
        """Find a file object corresponding to an include path."""
        for testdir in ('src', 'src/external/thread_mpi/include',
                'src/external/tng_io/include'):
            testpath = os.path.join(testdir, includedpath)
            if testpath in self._files:
                return self._files[testpath]

    def load_git_attributes(self):
        """Load git attribute information for files."""
        args = ['git', 'check-attr', '--stdin', 'filter']
        git_check_attr = subprocess.Popen(args, stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, cwd=self._source_root)
        filelist = '\n'.join(map(File.get_relpath, self._files.itervalues()))
        filters = git_check_attr.communicate(filelist)[0]
        for fileinfo in filters.splitlines():
            path, dummy, value = fileinfo.split(': ')
            fileobj = self._files.get(path)
            assert fileobj is not None
            fileobj.set_git_filter_attribute(value)

    def find_define_file_uses(self):
        """Find files that use defines from config.h."""
        # Executing git grep is substantially faster than using the define_re
        # directly on the contents of the file in Python.
        for define_file in self.get_checked_define_files():
            excluded_files = set([define_file])
            excluded_files.update(define_file.get_included_files(recursive=True))
            all_defines = define_file.get_declared_defines()
            args = ['git', 'grep', '-zwIF']
            for define in all_defines:
                args.extend(['-e', define])
            args.extend(['--', '*.cpp', '*.c', '*.cu', '*.h', '*.cuh'])
            define_re = r'\b(?:' + '|'.join(all_defines)+ r')\b'
            output = subprocess.check_output(args, cwd=self._source_root)
            for line in output.splitlines():
                (filename, text) = line.split('\0')
                fileobj = self._files.get(filename)
                if fileobj is not None and fileobj not in excluded_files:
                    defines = re.findall(define_re, text)
                    fileobj.add_used_defines(define_file, defines)

    def load_installed_file_list(self):
        """Load list of installed files from the build tree."""
        listpath = os.path.join(self._build_root, 'src', 'gromacs', 'installed-headers.txt')
        with open(listpath, 'r') as installedfp:
            for line in installedfp:
                path = line.strip()
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

    def load_cycle_suppression_list(self, filename):
        """Load a list of edges to suppress in cycles.

        These edges between modules, if present, will be marked in the
        corresponding ModuleDependency objects.
        """
        with open(filename, 'r') as fp:
            for line in fp:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                modulenames = ['module_' + x.strip() for x in line.split('->')]
                if len(modulenames) != 2:
                    self._reporter.input_error(
                            "invalid cycle suppression line: {0}".format(line))
                    continue
                firstmodule = self._modules.get(modulenames[0])
                secondmodule = self._modules.get(modulenames[1])
                if not firstmodule or not secondmodule:
                    self._reporter.input_error(
                            "unknown modules mentioned on cycle suppression line: {0}".format(line))
                    continue
                for dep in firstmodule.get_dependencies():
                    if dep.get_other_module() == secondmodule:
                        # TODO: Check that each suppression is actually part of
                        # a cycle.
                        dep.set_cycle_suppression()

    def get_object(self, docobj):
        """Get tree object for a Doxygen XML object."""
        if docobj is None:
            return None
        return self._docmap.get(docobj)

    def get_files(self):
        """Get iterable for all files in the source tree."""
        return self._files.itervalues()

    def get_modules(self):
        """Get iterable for all modules in the source tree."""
        return self._modules.itervalues()

    def get_classes(self):
        """Get iterable for all classes in the source tree."""
        return self._classes

    def get_members(self):
        """Get iterable for all members (in Doxygen terms) in the source tree."""
        return self._members

    def get_checked_define_files(self):
        """Get list of files that contain #define macros whose usage needs to
        be checked."""
        return (self._files['src/config.h'],
                self._files['src/gromacs/simd/simd.h'],
                self._files['src/gromacs/ewald/pme-simd.h'],
                self._files['src/gromacs/mdlib/nbnxn_simd.h'])
