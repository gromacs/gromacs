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

from gmxtree import DocType

class EdgeType(object):

    _names = ['undocumented', 'legacy', 'intramodule', 'public', 'library',
            'libimpl', 'pubimpl', 'test']

    def __init__(self, value):
        self._value = value

    def __str__(self):
        """Return string representation for the edge type (for debugging)."""
        return self._names[self._value]

    def __cmp__(self, other):
        """Order documentation types in the order of visibility."""
        return cmp(self._value, other._value)

# Tests depend on test
EdgeType.test = EdgeType(0)
# Implementation depends on public/library headers
EdgeType.pubimpl = EdgeType(1)
EdgeType.libimpl = EdgeType(2)
# Library header depends on other module
EdgeType.library = EdgeType(3)
# Public header depends on other module
EdgeType.public = EdgeType(4)
# Intramodule dependency
EdgeType.intramodule = EdgeType(5)
EdgeType.legacy = EdgeType(6)
# Invalid dependency
EdgeType.undocumented = EdgeType(7)

class Edge(object):

    """Graph edge between two node objects.

    Signifies an include dependency between the two nodes, and manages types
    associated with the dependencies.
    """

    def __init__(self, fromnode, tonode, edgetype):
        self._fromnode = fromnode
        self._tonode = tonode
        self._edgetype = edgetype

    def merge_edge(self, other):
        """Merge another edge into this one and choose an appropriate type.

        Updates the type of this edge based on the types of the merged edges.
        """
        self._edgetype = max(self._edgetype, other._edgetype)

    def format(self):
        """Format this edge for 'dot'."""
        if self._fromnode.is_file_node() and self._tonode.is_file_node():
            properties = ''
        elif self._edgetype == EdgeType.intramodule:
            properties = ''
        elif self._edgetype == EdgeType.test:
            properties = 'color=".33 .8 .8", style=dashed'
        elif self._edgetype == EdgeType.libimpl:
            properties = 'color=".66 .8 .8", style=dashed'
        elif self._edgetype == EdgeType.pubimpl:
            properties = 'color=black, style=dashed'
        elif self._edgetype == EdgeType.library:
            properties = 'color=".66 .8 .8"'
        elif self._edgetype == EdgeType.public:
            properties = 'color=black'
        elif self._edgetype == EdgeType.legacy:
            properties = 'color=grey75'
        else: # undocumented
            properties = 'color=red'
        return '{0} -> {1} [{2}]'.format(self._fromnode.get_nodename(),
                                         self._tonode.get_nodename(),
                                         properties)

class Node(object):
    def __init__(self, nodename, label, properties, is_file):
        self._nodename = nodename
        self._label = label
        self._properties = properties
        self._is_file = is_file
        self._children = []

    def add_child(self, child):
        self._children.append(child)

    def remove_child(self, child):
        self._children.remove(child)

    def clear_children(self):
        self._children = []

    def is_file_node(self):
        return self._is_file

    def get_nodename(self):
        return self._nodename

    def get_children(self, recursive=False):
        if recursive:
            result = list(self._children)
            for child in self._children:
                result.extend(child.get_children(recursive=True))
            return result
        else:
            return self._children

    def format(self):
        """Format this node for 'dot'."""
        result = ''
        if self._children:
            result += '    subgraph cluster_{0} {{\n' \
                          .format(self._nodename)
            result += '        label = "{0}"\n'.format(self._label)
            for child in self._children:
                result += child.format()
            result += '    }\n'
        else:
            properties = 'label="{0}"'.format(self._label)
            if self._properties:
                properties += ', ' + self._properties
            result += '    {0} [{1}]\n'.format(self._nodename, properties)
        return result


class Graph(object):
    def __init__(self, nodes, edges):
        self._nodes = set(nodes)
        self._edges = edges
        self._left_to_right = False
        self._concentrate = True

    def set_options(self, left_to_right=None, concentrate=None):
        if left_to_right != None:
            self._left_to_right = left_to_right
        if concentrate != None:
            self._concentrate = concentrate

    def merge_nodes(self, nodes, target):
        nodes = set(nodes)
        nodes.add(target)
        newedges = []
        edgesto = dict()
        edgesfrom = dict()
        for edge in self._edges:
            isfrom = (edge._fromnode in nodes)
            isto = (edge._tonode in nodes)
            if isfrom and isto:
                pass
            elif isfrom:
                if not edge._tonode in edgesfrom:
                    edgesfrom[edge._tonode] = \
                            Edge(target, edge._tonode, edge._edgetype)
                else:
                    edgesfrom[edge._tonode].merge_edge(edge)
            elif isto:
                if not edge._fromnode in edgesto:
                    edgesto[edge._fromnode] = \
                            Edge(edge._fromnode, target, edge._edgetype)
                else:
                    edgesto[edge._fromnode].merge_edge(edge)
            else:
                newedges.append(edge)
        newedges.extend(edgesfrom.values())
        newedges.extend(edgesto.values())
        self._edges = newedges

    def collapse_node(self, node):
        nodes = node.get_children(recursive=True)
        self.merge_nodes(nodes, node)
        node.clear_children()

    def write(self, outfile):
        outfile.write('digraph includedeps {\n')
        if self._left_to_right:
            outfile.write('    rankdir = LR\n')
        if self._concentrate:
            outfile.write('    concentrate = true\n')
        outfile.write('    node [fontname="FreeSans",fontsize=10,height=.2,'
                                 'shape=box]\n')
        for node in self._nodes:
            outfile.write(node.format())
        for edge in self._edges:
            outfile.write('    ' + edge.format() + '\n')
        outfile.write('}\n')

class GraphBuilder(object):
    def __init__(self, tree):
        self._tree = tree

    def _create_file_node(self, fileobj, filenodes):
        nodename = re.subn(r'[-./]', '_', fileobj.get_relpath())[0]
        properties = []
        style = []
        properties.append('URL="\\ref {0}"'.format(fileobj.get_name()))
        if not fileobj.get_module():
            style.append('bold')
            properties.append('color=red')
        if fileobj.is_test_file():
            style.append('filled')
            properties.append('fillcolor=".33 .2 1"')
        elif fileobj.is_source_file():
            style.append('filled')
            properties.append('fillcolor=grey75')
        elif fileobj.get_api_type() == DocType.public:
            style.append('filled')
            properties.append('fillcolor=".66 .2 1"')
        elif fileobj.get_api_type() == DocType.library:
            style.append('filled')
            properties.append('fillcolor=".66 .5 1"')
        if style:
            properties.append('style="{0}"'.format(','.join(style)))
        node = Node(nodename, fileobj.get_name(), ', '.join(properties), is_file=True)
        filenodes[fileobj] = node
        return node

    def _create_file_edge(self, fromfile, tofile, filenodes):
        intramodule = (fromfile.get_module() == tofile.get_module())
        is_legacy = not tofile.get_module().is_documented()
        if fromfile.get_module() == tofile.get_module():
            edgetype = EdgeType.intramodule
        elif tofile.get_api_type() == DocType.internal:
            if is_legacy:
                edgetype = EdgeType.legacy
            else:
                edgetype = EdgeType.undocumented
        elif fromfile.is_test_file():
            edgetype = EdgeType.test
        elif tofile.is_test_file():
            edgetype = EdgeType.undocumented
        elif fromfile.is_source_file() or \
                (fromfile.get_api_type() <= DocType.internal and \
                not fromfile.is_installed()):
            if tofile.get_api_type() == DocType.public:
                edgetype = EdgeType.pubimpl
            elif tofile.get_api_type() == DocType.library:
                edgetype = EdgeType.libimpl
            elif is_legacy or not tofile.is_documented():
                edgetype = EdgeType.legacy
            else:
                raise ValueError('Unknown edge type between {0} and {1}'
                        .format(fromfile.path, tofile.path))
        elif fromfile.get_api_type() == DocType.library:
            edgetype = EdgeType.library
        elif fromfile.get_api_type() == DocType.public or fromfile.is_installed():
            if tofile.get_api_type() == DocType.public or \
                    tofile.get_documentation_type() == DocType.public or \
                    (tofile.is_installed() and not tofile.is_documented()):
                edgetype = EdgeType.public
            else:
                edgetype = EdgeType.undocumented
        else:
            raise ValueError('Unknown edge type between {0} and {1}'
                    .format(fromfile.path, tofile.path))
        return Edge(filenodes[fromfile], filenodes[tofile], edgetype)

    def _create_file_edges(self, filenodes):
        edges = []
        for fileobj in filenodes.iterkeys():
            for includedfile in fileobj.get_includes():
                otherfile = includedfile.get_file()
                if otherfile and otherfile in filenodes:
                    edge = self._create_file_edge(fileobj, otherfile, filenodes)
                    edges.append(edge)
        return edges

    def _create_module_node(self, module, filenodes):
        properties = 'shape=ellipse, URL="\\ref module_{0}"'.format(module.get_name())
        if not module.is_documented():
            properties += ', style=filled, fillcolor=grey75'
        elif module.get_group() == 'analysismodules':
            properties += ', style=filled, fillcolor="0 .2 1"'
        elif module.get_group() == 'utilitymodules':
            properties += ', style=filled, fillcolor=".08 .2 1"'
        elif module.get_group() == 'mdrun':
            properties += ', style=filled, fillcolor=".75 .2 1"'
        rootdir = module.get_root_dir()
        if rootdir.has_installed_files():
            properties += ', color=".66 .5 1", penwidth=3'
        nodename = 'module_' + re.subn(r'[-./]', '_', rootdir.get_name())[0]
        label = module.get_name()[7:]
        node = Node(nodename, label, properties, is_file=False)
        for childfile in module.get_files():
            node.add_child(self._create_file_node(childfile, filenodes))
        return node

    def create_modules_graph(self):
        filenodes = dict()
        nodes = []
        modulenodes = []
        libgromacsnode = Node('libgromacs', 'libgromacs', '', is_file=False)
        nodes.append(libgromacsnode)
        for moduleobj in self._tree.get_modules():
            node = self._create_module_node(moduleobj, filenodes)
            if moduleobj.get_root_dir().get_relpath().startswith('src/gromacs'):
                libgromacsnode.add_child(node)
            else:
                nodes.append(node)
            modulenodes.append(node)
        edges = self._create_file_edges(filenodes)
        graph = Graph(nodes, edges)
        for node in modulenodes:
            graph.collapse_node(node)
        graph.set_options(concentrate=False)
        return graph

    def create_module_file_graph(self, module):
        filenodes = dict()
        nodes = []
        for fileobj in module.get_files():
            nodes.append(self._create_file_node(fileobj, filenodes))
        edges = self._create_file_edges(filenodes)
        graph = Graph(nodes, edges)
        graph.set_options(left_to_right=True)
        return graph

def main():
    """Run the graph generation script."""
    import sys

    from optparse import OptionParser

    from gmxtree import GromacsTree
    from reporter import Reporter

    parser = OptionParser()
    parser.add_option('-S', '--source-root',
                      help='Source tree root directory')
    parser.add_option('-B', '--build-root',
                      help='Build tree root directory')
    parser.add_option('--installed',
                      help='Read list of installed files from given file')
    parser.add_option('-o', '--outdir', default='.',
                      help='Specify output directory for graphs')
    parser.add_option('-q', '--quiet', action='store_true',
                      help='Do not write status messages')
    options, args = parser.parse_args()

    installedlist = []
    if options.installed:
        with open(options.installed, 'r') as outfile:
            for line in outfile:
                installedlist.append(line.strip())

    reporter = Reporter(quiet=True)

    if not options.quiet:
        sys.stderr.write('Scanning source tree...\n')
    tree = GromacsTree(options.source_root, options.build_root, reporter)
    tree.set_installed_file_list(installedlist)
    if not options.quiet:
        sys.stderr.write('Reading source files...\n')
    tree.scan_files()
    if not options.quiet:
        sys.stderr.write('Reading Doxygen XML files...\n')
    tree.load_xml(only_files=True)

    if not options.quiet:
        sys.stderr.write('Writing graphs...\n')
    graphbuilder = GraphBuilder(tree)

    filename = os.path.join(options.outdir, 'module-deps.dot')
    graph = graphbuilder.create_modules_graph()
    with open(filename, 'w') as outfile:
        graph.write(outfile)

    # Skip some modules that are too big to make any sense
    skippedmodules = ('legacyheaders', 'gmxlib', 'mdlib', 'gmxana', 'gmxpreprocess')
    for module in tree.get_modules():
        if not module.get_name()[7:] in skippedmodules:
            filename = '{0}-deps.dot'.format(module.get_name())
            filename = os.path.join(options.outdir, filename)
            graph = graphbuilder.create_module_file_graph(module)
            with open(filename, 'w') as outfile:
                graph.write(outfile)

if __name__ == '__main__':
    main()
