#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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

"""Check and generate include dependency graphs for Gromacs.

This script can do a few things related to include file dependencies:
 - Check that there are no broken dependencies between installed headers.
 - Check that documentated usage of a header matches its installation status
   and usage from other modules.
 - Generate two types of include dependency graphs: per-file or per-module
   (where module is equivalent to a subdirectory).
It is intended to be run on a subset of files under the src/ directory.
Output format for the graphs is suitable for processing with 'dot'.

FILE GRAPHS

The legend for per-file graph nodes:
    gray:          source files
    light blue:    public headers
    dark blue:     library headers
    no background: other files

MODULE GRAPHS

Module graph will contain one node for each top-level subdirectory under src/,
except that the src/gromacs/ directory will be expanded one level further.
Legacy modules have gray background.

The legend for per-module graph links (a link with a certain color indicates
that types above it in the list are not present):
    red:          invalid dependency
    grey:         legacy dependency (dependency on undocumented file, or to
                  legacy directories)
    solid blue:   public header depends on the other module
    solid black:  library header depends on the other module
    dashed blue:  source file depends on a library header in the other module
    dashed black: source file depends on a public header in the other module
    dashed green: test file depends on the other module
"""

import os.path
import re
import sys

from optparse import OptionParser

class ErrorReporter(object):
    def input_warning(self, file_path, msg):
        sys.stderr.write('warning: {0}: {1}\n'.format(file_path, msg))

    def error(self, file_path, msg):
        sys.stderr.write('error: {0}: {1}\n'.format(file_path, msg))

class Link(object):

    """Link between two node objects.

    Signifies an include dependency between the two nodes, and manages types
    associated with the dependencies.
    """

    _priorities = {
            'undocumented': 1,
            'legacy': 2,
            'intramodule': 3,
            'public': 4,
            'library': 5,
            'libimpl': 6,
            'pubimpl': 7,
            'test': 8}

    def __init__(self, fromnode, tonode, link_type):
        self.fromnode = fromnode
        self.tonode = tonode
        self.link_type = link_type
        if link_type not in Link._priorities:
            raise ValueError('Unknown link type {0}'.format(link_type))

    def merge_link(self, other):
        """Merge another link into this one and choose an appropriate type.

        Updates the type of this link based on the types of the merged links.
        """
        if Link._priorities[other.link_type] < Link._priorities[self.link_type]:
            self.link_type = other.link_type

    def format(self):
        """Format this link for 'dot'."""
        if self.fromnode.is_file_node() and self.tonode.is_file_node():
            properties = ''
        elif self.link_type == 'intramodule':
            properties = ''
        elif self.link_type == 'test':
            properties = 'color=".33 .8 .8", style=dashed'
        elif self.link_type == 'libimpl':
            properties = 'color=".66 .8 .8", style=dashed'
        elif self.link_type == 'pubimpl':
            properties = 'color=black, style=dashed'
        elif self.link_type == 'library':
            properties = 'color=".66 .8 .8"'
        elif self.link_type == 'public':
            properties = 'color=black'
        elif self.link_type == 'legacy':
            properties = 'color=grey75'
        else: # undocumented
            properties = 'color=red'
        return '{0} -> {1} [{2}]'.format(self.fromnode.nodename,
                                         self.tonode.nodename,
                                         properties)

class Node(object):
    def __init__(self, nodename, label, properties, is_file):
        self.nodename = nodename
        self.label = label
        self._properties = properties
        self._is_file = is_file
        self.children = []
        self.root = False

    def set_root(self):
        self.root = True

    def add_child(self, child):
        self.children.append(child)

    def remove_child(self, child):
        self.children.remove(child)

    def clear_children(self):
        self.children = []

    def is_file_node(self):
        return self._is_file

    def get_children(self, recursive=False):
        if recursive:
            result = list(self.children)
            for child in self.children:
                result.extend(child.get_children(recursive=True))
            return result
        else:
            return self.children

    def format(self):
        """Format this node for 'dot'."""
        result = ''
        if self.children:
            if not self.root:
                result += '    subgraph cluster_{0} {{\n' \
                              .format(self.nodename)
                result += '        label = "{0}"\n'.format(self.label)
            for child in self.children:
                result += child.format()
            if not self.root:
                result += '    }\n'
        else:
            properties = 'label="{0}"'.format(self.label)
            if self._properties:
                properties += ', ' + self._properties
            result += '    {0} [{1}]\n'.format(self.nodename, properties)
        return result


class Graph(object):
    def __init__(self, nodes, links):
        self.nodes = set(nodes)
        self.links = links
        self.left_to_right = False
        self.concentrate = True

    def set_options(self, left_to_right=None, concentrate=None):
        if left_to_right != None:
            self.left_to_right = left_to_right
        if concentrate != None:
            self.concentrate = concentrate

    def prune_links(self):
        nodes = set()
        for node in self.nodes:
            nodes.update(node.get_children(recursive=True))
        newlinks = []
        for link in self.links:
            if link.fromnode in nodes and link.tonode in nodes:
                newlinks.append(link)
        self.links = newlinks

    def merge_nodes(self, nodes, target):
        nodes = set(nodes)
        nodes.add(target)
        newlinks = []
        linksto = dict()
        linksfrom = dict()
        for link in self.links:
            isfrom = (link.fromnode in nodes)
            isto = (link.tonode in nodes)
            if isfrom and isto:
                pass
            elif isfrom:
                if not link.tonode in linksfrom:
                    linksfrom[link.tonode] = \
                            Link(target, link.tonode, link.link_type)
                else:
                    linksfrom[link.tonode].merge_link(link)
            elif isto:
                if not link.fromnode in linksto:
                    linksto[link.fromnode] = \
                            Link(link.fromnode, target, link.link_type)
                else:
                    linksto[link.fromnode].merge_link(link)
            else:
                newlinks.append(link)
        newlinks.extend(linksfrom.values())
        newlinks.extend(linksto.values())
        self.links = newlinks

    def collapse_node(self, node):
        nodes = node.get_children(recursive=True)
        self.merge_nodes(nodes, node)
        node.clear_children()

    def write(self, outfile):
        outfile.write('digraph includedeps {\n')
        if self.left_to_right:
            outfile.write('    rankdir = LR\n')
        if self.concentrate:
            outfile.write('    concentrate = true\n')
        outfile.write('    node [fontname="FreeSans",fontsize=10,height=.2,'
                                 'shape=box]\n')
        for link in self.links:
            outfile.write('    ' + link.format() + '\n')
        for node in self.nodes:
            outfile.write(node.format())
        outfile.write('}\n')


def find_include_file(filename, includedirs):
    """Find full path to filename, looking in a set of directories."""
    for includedir in includedirs:
        fullpath = os.path.abspath(os.path.join(includedir, filename))
        if os.path.exists(fullpath):
            return fullpath
    return None


class IncludedFile(object):
    def __init__(self, included_file, included_path, is_relative, is_system):
        self._included_file = included_file
        self._included_path = included_path
        #self._used_include_path = used_include_path
        self._is_relative = is_relative
        self._is_system = is_system


class File(object):
    def __init__(self, path, module):
        self.path = path
        self.name = os.path.basename(path)
        self.module = module
        if module.name == 'tests':
            self.type = 'test'
        elif re.search(r'\.c(pp|u)?$', self.name) != None:
            self.type = 'source'
        else:
            self.type = 'header'
        self.doctype = 'none'
        #headername = re.sub(r'\.cpp$', '.h', self.name)
        #implheadername = re.sub(r'\.cpp$', '-impl.h', self.name)
        self._included = []
        self.installed = False

    def is_documented(self):
        return self.doctype != 'none'

    def is_installed(self):
        return self.installed

    def set_installed(self, reporter):
        if self.type != 'header':
            reporter.input_warning(self.path,
                    'installing {0} file'.format(self.type))
            return
        self.installed = True

    def get_included_files(self):
        return self._included

    def scan_include_file(self, line, allfiles, selfdir, includedirs,
            ignorelist, reporter):
        """Process #include directive during scan().

        Searches for the included file in given directories, does some checks,
        and adds the dependency link to the other file if applicable.
        """
        fullpath = None
        includedpath = None
        includedfile = None
        is_system = False
        is_relative = False
        match = re.match(r'#include *<([^>]*)>', line)
        if match:
            includedpath = match.group(1)
            is_system = True
            fullpath = find_include_file(includedpath, includedirs)
        else:
            match = re.match(r'#include *"([^"]*)"', line)
            if match:
                includedpath = match.group(1)
                fullpath = os.path.abspath(os.path.join(selfdir, includedpath))
                #if os.path.abspath(fullpath) in ignorelist:
                #    return
                if os.path.exists(fullpath):
                    is_relative = True
                else:
                    fullpath = find_include_file(includedpath, includedirs)
                    if not fullpath:
                        if not includedpath in ('corewrap.h', 'tmpi_config.h'):
                            reporter.input_warning(self.path,
                                    'included file "{0}" not found'
                                        .format(includedpath))
        if not includedpath:
            reporter.input_warning(self.path, 'line "{0}" could not be parsed'
                    .format(line))
        else:
            if fullpath and fullpath in allfiles:
                includedfile = allfiles[fullpath]
            #elif not dep in ignorelist:
            #    depfile = File(dep, None)
            #    files[dep] = depfile
            #    file.add_dependency(depfile)
            #    extrafiles.append(dep)
            self._included.append(IncludedFile(includedfile, includedpath,
                    is_relative, is_system))

    def scan(self, filename, allfiles, includedirs, ignorelist, reporter):
        selfdir = os.path.dirname(filename)
        infileblock = False
        foundfileblock = False
        self.docmodule = None
        with open(filename, 'r') as scanfile:
            for line in scanfile:
                if line.startswith('#include'):
                    self.scan_include_file(line, allfiles, selfdir,
                            includedirs, ignorelist, reporter)
                    continue
                if not foundfileblock:
                    if infileblock:
                        if r'*/' in line:
                            infileblock = False
                            foundfileblock = True
                            continue
                        if self.type == 'implheader':
                            if line.startswith(r' * \inpublicapi'):
                                self.type = 'publicheader'
                            elif line.startswith(r' * \inlibraryapi'):
                                self.type = 'libheader'
                        match = re.match(r' \* \\ingroup module_([a-z_]*)', line)
                        if match:
                            if self.docmodule:
                                reporter.error(self.path,
                                        'file documented in multiple modules')
                            self.docmodule = match.group(1)
                    else:
                        match = re.match(r'/\*! *(\\[a-z]*internal)? *\\file', line)
                        if match:
                            docspec = match.group(1)
                            if not docspec:
                                self.doctype = 'public'
                            elif docspec == r'\libinternal':
                                self.doctype = 'library'
                            elif docspec == r'\internal':
                                self.doctype = 'implementation'
                            else:
                                reporter.input_warning(self.path,
                                        'unknown specifier "{0}"'.format(docspec))
                                self.doctype = 'unknown'
                            infileblock = True
                            if self.type == 'header':
                                # Default type if no other found
                                self.type = 'implheader'


class Module(object):
    def __init__(self, name, parent = None):
        self.parent = parent
        self.name = name
        if parent:
            self.fullname = parent.fullname + '_' + name
        else:
            self.fullname = 'module'
        self.files = []
        self.children = dict()
        self.is_top_level = (not parent or parent.name in ('', 'gromacs'))

    def get_parent(self):
        return self.parent

    def is_child(self, module):
        parent = module.parent
        while parent:
            if parent == self:
                return True
            parent = parent.parent
        return False

    def get_top_level_module(self):
        if self.is_top_level or not self.parent:
            return self
        return self.parent.get_top_level_module()

    def add_nested_file(self, modules, path):
        if len(modules) == 1:
            newfile = File(path, self)
            self.files.append(newfile)
        else:
            if not modules[0] in self.children:
                module = Module(modules[0], self)
                self.children[modules[0]] = module
            else:
                module = self.children[modules[0]]
            newfile = module.add_nested_file(modules[1:], path)
        return newfile


class Dependencies(object):
    def __init__(self, rootdir, includedirs, installedfiles):
        self.files = dict()
        self.root = Module("")
        self.rootpath = []
        for root in rootdir:
            self.rootpath.append(os.path.abspath(root))
        if includedirs:
            self.includedirs = self.rootpath + includedirs
        else:
            self.includedirs = self.rootpath
        self.installedfiles = installedfiles

    def add_file(self, filename, reporter):
        fullpath = os.path.abspath(filename)
        for root in self.rootpath:
            if fullpath.startswith(root):
                relpath = fullpath[len(root)+1:]
                break
        else:
            reporter.input_warning(filename,
                    'input file not under root path, skipped')
            return
        modules = relpath.split(os.sep)
        newfile = self.root.add_nested_file(modules, relpath)
        if fullpath in self.installedfiles:
            newfile.set_installed(reporter)
        self.files[os.path.abspath(filename)] = newfile

    def scan_files(self, ignorelist, reporter):
        for (filename, scanfile) in self.files.iteritems():
            scanfile.scan(filename, self.files, self.includedirs, ignorelist,
                    reporter)

    def get_toplevel_modules(self):
        result = []
        for module in self.root.children.itervalues():
            if module.name == 'gromacs':
                result.extend(module.children.itervalues())
            else:
                result.append(module)
        return result


def _is_legacy_module(module):
    if module.name in ('legacyheaders', 'gmxlib', 'mdlib', 'gmxana', 'gmxpreprocess'):
        return True
    if module.get_parent():
        return _is_legacy_module(module.get_parent())
    return False


class IncludeFileChecker(object):
    def __init__(self, deps, options):
        self._deps = deps
        self._options = options

    def _check_file(self, checkfile, reporter):
        if not self._options.check_doc:
            return
        if not checkfile.is_documented():
            if self._options.warn_undoc:
                is_legacy = _is_legacy_module(checkfile.module)
                is_external = checkfile.module.name in ('gmx_lapack', 'gmx_blas', 'thread_mpi')
                if not is_legacy and not is_external:
                    reporter.error(checkfile.path, 'file not documented')
        elif checkfile.doctype == 'implementation' and \
                checkfile.type in ('publicheader', 'libheader'):
            reporter.error(checkfile.path,
                    'file documentation visibility incorrect')
        elif checkfile.doctype == 'library' and checkfile.type == 'publicheader':
            reporter.error(checkfile.path,
                    'file documentation visibility incorrect')
        elif checkfile.installed and checkfile.doctype not in ('public', 'unknown'):
            reporter.error(checkfile.path,
                    'installed header has no public documentation')
        elif not checkfile.installed and checkfile.doctype == 'public':
            reporter.error(checkfile.path,
                    'non-installed file has public documentation')
        selfmodfullname = checkfile.module.fullname
        docmodule = checkfile.docmodule
        if docmodule and \
                not selfmodfullname.startswith('module_' + docmodule) and \
                not selfmodfullname.startswith('module_gromacs_' + docmodule) and \
                not checkfile.name == docmodule + '.h':
            reporter.error(checkfile.path,
                    'file documented in incorrect module "{0}"'
                        .format(docmodule))

    def _check_included_file(self, checkfile, includedfile, reporter):
        otherfile = includedfile._included_file
        if includedfile._is_system:
            # TODO: This doesn't report errors with files not listed in
            # the input files, although those could be included.
            # That would produce a massive amount of errors for <config.h>.
            if otherfile:
                reporter.error(checkfile.path,
                        'local file included as <{0}>'
                            .format(includedfile._included_path))
        elif not includedfile._is_relative and checkfile.installed:
            if not includedfile._included_path == 'gmx_header_config_gen.h':
                reporter.error(checkfile.path,
                        'installed header includes "{0}", '
                        'which is not found using relative path'
                            .format(includedfile._included_path))
        if not otherfile:
            return
        if checkfile.installed and not otherfile.installed:
            reporter.error(checkfile.path,
                    'installed header includes '
                    'non-installed header "{0}"'
                        .format(includedfile._included_path))
        if not otherfile.is_documented():
            return
        if not self._options.check_doc:
            return
        intramodule = \
                (checkfile.module.get_top_level_module() == \
                 otherfile.module.get_top_level_module())
        if otherfile.type not in ('publicheader', 'libheader'):
            if not intramodule and not _is_legacy_module(otherfile.module):
                reporter.error(checkfile.path,
                        'included file "{0}" is missing API definition'
                            .format(otherfile.path))
        elif checkfile.type == 'publicheader':
            if not otherfile.type == 'publicheader' and not otherfile.doctype == 'public':
                reporter.error(checkfile.path,
                        'public API file includes non-public header "{0}"'
                            .format(otherfile.path))

    def check_all(self, reporter):
        for checkfile in sorted(self._deps.files.values()):
            self._check_file(checkfile, reporter)
            for includedfile in checkfile.get_included_files():
                self._check_included_file(checkfile, includedfile, reporter)


class GraphBuilder(object):
    def __init__(self, deps):
        self._deps = deps

    def _create_file_node(self, fileobj, filenodes):
        nodename = re.subn(r'[-./]', '_', fileobj.path)[0]
        properties = []
        style = []
        properties.append('URL="\\ref {0}"'.format(fileobj.name))
        if not fileobj.module:
            style.append('bold')
            properties.append('color=red')
        if fileobj.type == 'source':
            style.append('filled')
            properties.append('fillcolor=grey75')
        elif fileobj.type == 'publicheader':
            style.append('filled')
            properties.append('fillcolor=".66 .2 1"')
        elif fileobj.type == 'libheader':
            style.append('filled')
            properties.append('fillcolor=".66 .5 1"')
        if style:
            properties.append('style="{0}"'.format(','.join(style)))
        node = Node(nodename, fileobj.name, ', '.join(properties), is_file=True)
        filenodes[fileobj] = node
        return node

    def _create_file_edge(self, fromfile, tofile, filenodes):
        intramodule = \
                (fromfile.module.get_top_level_module() == \
                 tofile.module.get_top_level_module())
        is_legacy = _is_legacy_module(tofile.module)
        if tofile.type not in ('publicheader', 'libheader', 'header'):
            if intramodule:
                link_type = 'intramodule'
            elif is_legacy:
                link_type = 'legacy'
            else:
                link_type = 'undocumented'
        elif fromfile.type == 'test':
            link_type = 'test'
        elif fromfile.type in ('source', 'header', 'implheader') and \
                not fromfile.is_installed():
            if intramodule:
                link_type = 'intramodule'
            elif tofile.type == 'publicheader':
                link_type = 'pubimpl'
            elif tofile.type == 'libheader':
                link_type = 'libimpl'
            elif is_legacy:
                link_type = 'legacy'
            elif not tofile.is_documented():
                link_type = 'legacy'
            else:
                raise ValueError('Unknown link type between {0} and {1}'
                        .format(fromfile.path, tofile.path))
        elif fromfile.type == 'libheader':
            link_type = 'library'
        elif fromfile.type == 'publicheader' or fromfile.is_installed():
            if tofile.type == 'publicheader' or tofile.doctype == 'public' or \
                    (tofile.is_installed() and not tofile.is_documented()):
                link_type = 'public'
            else:
                link_type = 'undocumented'
        else:
            raise ValueError('Unknown link type between {0} and {1}'
                    .format(fromfile.path, tofile.path))
        return Link(filenodes[fromfile], filenodes[tofile], link_type)

    def _create_file_edges(self, fileobj, filenodes):
        links = []
        if fileobj in filenodes:
            for includedfile in fileobj.get_included_files():
                otherfile = includedfile._included_file
                if otherfile and otherfile in filenodes:
                    link = self._create_file_edge(fileobj, otherfile, filenodes)
                    links.append(link)
        return links

    def create_module_node(self, module, filenodes):
        properties = 'shape=ellipse, URL="\\ref module_{0}"'.format(module.name)
        if _is_legacy_module(module):
            properties += 'style=filled, fillcolor=grey75'
        node = Node(module.fullname, module.name, properties, is_file=False)
        for childfile in module.files:
            node.add_child(self._create_file_node(childfile, filenodes))
        for childmodule in module.children.itervalues():
            node.add_child(self.create_module_node(childmodule, filenodes))
        return node

    def create_file_graph(self):
        filenodes = dict()
        rootnode = self.create_module_node(self._deps.root, filenodes)
        rootnode.set_root()
        links = []
        for scanfile in self._deps.files.itervalues():
            links.extend(self._create_file_edges(scanfile, filenodes))
        graph = Graph([rootnode], links)
        return graph

    def create_modules_graph(self):
        filenodes = dict()
        rootnode = self.create_module_node(self._deps.root, filenodes)
        rootnode.set_root()
        links = []
        for scanfile in self._deps.files.itervalues():
            links.extend(self._create_file_edges(scanfile, filenodes))
        graph = Graph([rootnode], links)
        for node in rootnode.get_children():
            if node.label == 'gromacs':
                module_nodes = []
                header_nodes = []
                for child in node.get_children():
                    if child.is_file_node():
                        header_nodes.append(child)
                    else:
                        graph.collapse_node(child)
                        module_nodes.append(child)
                for header in header_nodes:
                    for module in module_nodes:
                        if header.nodename.startswith(module.nodename[7:]):
                            # graph.merge_nodes([header], module)
                            node.remove_child(header)
                            break
            else:
                graph.collapse_node(node)
        graph.set_options(concentrate=False)
        graph.prune_links()
        return graph

    def create_module_file_graph(self, module):
        filenodes = dict()
        rootnode = self.create_module_node(module, filenodes)
        rootnode.set_root()
        links = []
        for scanfile in self._deps.files.itervalues():
            links.extend(self._create_file_edges(scanfile, filenodes))
        graph = Graph([rootnode], links)
        graph.prune_links()
        return graph


def print_module_graph(outfile, graphbuilder, options):
    graph = graphbuilder.create_modules_graph()
    graph.write(outfile)

def print_file_graph(outfile, graphbuilder, options):
    graph = graphbuilder.create_file_graph()
    graph.set_options(left_to_right=options.left_to_right)
    graph.write(outfile)
    #if options.source_at_top:
    #    sourcenodes = []
    #    for file in deps.files.itervalues():
    #        if file.sourcefile:
    #            sourcenodes.append(file.nodename)
    #    if sourcenodes:
    #        outfile.write('    { rank = min; ' + '; '.join(sourcenodes) + '}\n')
    #if options.with_external and options.external_at_bottom:
    #    extnodes = []
    #    for file in deps.files.itervalues():
    #        if not file.module:
    #            extnodes.append(file.nodename)
    #    if extnodes:
    #        outfile.write('    { rank = max; ' + '; '.join(extnodes) + '}\n')

def print_module_file_graph(outfile, graphbuilder, module, options):
    graph = graphbuilder.create_module_file_graph(module)
    graph.set_options(left_to_right=options.left_to_right)
    graph.write(outfile)

def main():
    parser = OptionParser()
    parser.add_option('-f', '--files',
                      help='Read list of input files from given file')
    parser.add_option('--installed',
                      help='Read list of installed files from given file')
    parser.add_option('-R', '--rootdir', action='append',
                      help='Remove this prefix from all files')
    parser.add_option('-I', '--includedir', action='append',
                      help='Specify additional directories to search for '
                           'include files')
    parser.add_option('-o', '--outdir', default='.',
                      help='Specify output directory for graphs')
    #parser.add_option('--source-at-top', action='store_true',
    #                  help='Force source files at the top of the graph')
    #parser.add_option('--with-external', action='store_true',
    #                  help='Include external dependencies in the graph')
    #parser.add_option('--external-at-bottom', action='store_true',
    #                  help='Force external dependencies files at the bottom '
    #                       'of the graph')
    parser.add_option('--check', action='store_true',
                      help='Check for problems in include file dependencies')
    parser.add_option('--check-doc', action='store_true',
                      help='Check for problems in Doxygen documentation')
    parser.add_option('--warn-undoc', action='store_true',
                      help='Warn for files that do not have Doxygen documentation')
    parser.add_option('--left-to-right', action='store_true',
                      help='Lay out from left to right')
    parser.add_option('--file-graph',
                      help='Write graph for individual files')
    parser.add_option('--module-graph',
                      help='Write graph for modules')
    parser.add_option('--module-file-graphs', action='store_true',
                      help='Write file graphs for each module')
    options, args = parser.parse_args()

    if not options.file_graph and not options.module_graph and \
            not options.module_file_graphs:
        options.check = True

    # Constructs lists of files
    filelist = []
    ignorelist = []
    installedlist = []
    if options.files:
        with open(options.files, 'r') as outfile:
            for line in outfile:
                if line.startswith('!'):
                    ignorelist.append(os.path.abspath(line[1:].strip()))
                else:
                    filelist.append(line.strip())
    filelist.extend(args)
    if options.installed:
        with open(options.installed, 'r') as outfile:
            for line in outfile:
                installedlist.append(line.strip())

    # Creates objects for all files and modules
    reporter = ErrorReporter()
    deps = Dependencies(options.rootdir, options.includedir, installedlist)
    for filename in filelist:
        deps.add_file(filename, reporter)

    deps.scan_files(ignorelist, reporter)

    if options.check or options.check_doc:
        checker = IncludeFileChecker(deps, options)
        checker.check_all(reporter)

    #if options.with_external:
    #    for filename in extrafiles:
    #        file = files[filename]
    #        if os.path.exists(filename):
    #            with open(filename, 'r') as outfile:
    #                for line in outfile:
    #                    if not file.api:
    #                        if line.startswith(' * \inpublicapi'):
    #                            file.api = "public"
    #                        elif line.startswith(' * \inlibraryapi'):
    #                            file.api = "library"

    # Prints out the graph
    graphbuilder = GraphBuilder(deps)
    if options.module_graph:
        graphpath = os.path.join(options.outdir, options.module_graph)
        with open(graphpath, 'w') as outfile:
            print_module_graph(outfile, graphbuilder, options)
    if options.file_graph:
        graphpath = os.path.join(options.outdir, options.file_graph)
        with open(graphpath, 'w') as outfile:
            print_file_graph(outfile, graphbuilder, options)
    if options.module_file_graphs:
        options.left_to_right = True
        for module in deps.get_toplevel_modules():
            if not _is_legacy_module(module):
                filename = 'module_{0}-deps.dot'.format(module.name)
                filename = os.path.join(options.outdir, filename)
                with open(filename, 'w') as outfile:
                    print_module_file_graph(outfile, graphbuilder, module, options)

main()
