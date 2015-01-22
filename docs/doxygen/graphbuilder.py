#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

"""Generate include dependency graphs.

This script generates include dependency graphs from the GROMACS source tree.
One graph is generated to show inter-module dependencies, and separate graphs
for each module to show file-level dependencies within the module.

Output format for the graphs is suitable for processing with 'dot' in graphviz.

The graphs are built from the source tree representation constructed in
gmxtree.py.

Classes Graph, Node, Edge, and EdgeType provide a relatively general
implementation for constructing 'dot' graphs.  GraphBuilder is used to
create Graph instances from a gmxtree.GromacsTree object; the actual graph
objects will not contain any references to the gmxtree objects.

When run in script mode, the GromacsTree object is first constructed, and then
GraphBuilder is used to construct the necessary graphs, which are then written
out.

The produced graphs are documented in doxygen.md.
"""

import os.path
import re

from gmxtree import DocType

class EdgeType(object):

    """Enumeration type for edge types in include dependency graphs."""

    # Mapping to string representation for the internal integer values
    _names = ['test', 'pubimpl', 'libimpl', 'library', 'public',
            'intramodule', 'legacy', 'undocumented']

    def __init__(self, value):
        """Initialize a EdgeType instance.

        EdgeType.{test,pubimpl,...,undocumented} should be used outside the
        class instead of calling the constructor.
        """
        self._value = value

    def __str__(self):
        """Return string representation for the edge type (for debugging)."""
        return self._names[self._value]

    def __cmp__(self, other):
        """Order edge types in the order of increasing coupling."""
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
EdgeType.cyclic = EdgeType(7)
# Invalid dependency
EdgeType.undocumented = EdgeType(8)

class Edge(object):

    """Graph edge between two Node objects in 'dot' graph.

    Signifies an include dependency between the two nodes, and manages types
    associated with the dependencies.
    """

    def __init__(self, fromnode, tonode, edgetype):
        """Create edge between given Nodes with given type."""
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
        # If you change these styles, update also the legend in modulegraph.md
        if self._fromnode.is_file_node() and self._tonode.is_file_node():
            properties = ''
        elif self._edgetype == EdgeType.intramodule:
            properties = ''
        elif self._edgetype == EdgeType.test:
            # TODO: Consider if only some test edges should be made non-constraints
            properties = 'color=".33 .8 .8", style=dashed, constraint=no'
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
        elif self._edgetype == EdgeType.cyclic:
            properties = 'color=red, constraint=no'
        else: # undocumented
            properties = 'color=red'
        return '{0} -> {1} [{2}]'.format(self._fromnode.get_nodename(),
                                         self._tonode.get_nodename(),
                                         properties)

class Node(object):

    """Node in 'dot' graph."""

    def __init__(self, nodename, label, style=None, properties=None, is_file=False):
        """Create node with given attributes.

        is_file does not affect the appearance of the node, but is used for
        formatting edges between two files differently from other edges.
        style and properties should be iterables with graphviz attributes for
        the node.

        Node can have child nodes.  Such nodes are rendered as cluster
        subgraphs for 'dot'.
        """
        self._nodename = nodename
        self._label = label
        if style:
            self._style = ','.join(style)
        else:
            self._style = None
        if properties:
            self._properties = ', '.join(properties)
        else:
            self._properties = None
        self._is_file = is_file
        self._children = []

    def add_child(self, child):
        """Add a child node."""
        self._children.append(child)

    def clear_children(self):
        """Remove all children from the node."""
        self._children = []

    def is_file_node(self):
        """Return True if the node was created with is_file=True."""
        return self._is_file

    def get_nodename(self):
        """Get internal name of the node in 'dot'."""
        return self._nodename

    def get_children(self, recursive=False):
        """Get list of child nodes."""
        if recursive:
            result = list(self._children)
            for child in self._children:
                result.extend(child.get_children(recursive=True))
            return result
        else:
            return self._children

    def format(self):
        """Format this node for 'dot'."""
        # TODO: Take indent as a parameter to make output marginally nicer.
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
            if self._style:
                properties += ', style="{0}"'.format(self._style)
            result += '    {0} [{1}]\n'.format(self._nodename, properties)
        return result


class Graph(object):

    """Graph for 'dot'."""

    def __init__(self, nodes, edges):
        """Create graph with given nodes and edges."""
        self._nodes = set(nodes)
        self._edges = edges
        self._left_to_right = False
        self._concentrate = True

    def set_options(self, left_to_right=None, concentrate=None):
        """Set output options for the graph."""
        if left_to_right != None:
            self._left_to_right = left_to_right
        if concentrate != None:
            self._concentrate = concentrate

    def merge_nodes(self, nodes, target):
        """Merge a set of nodes into a single node.

        All nodes from the list nodes are merged into the target node.
        All edges to or from the merged nodes are rerouted to/from target
        instead.  Duplicate edges are not created.  Instead, if an edge already
        exists, the edge types are merged.  All nodes from the list nodes are
        removed from the graph after the merge is done.
        """
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
        """Merge all children of a node into the node.

        All child nodes are removed after the merge is done.
        """
        nodes = node.get_children(recursive=True)
        self.merge_nodes(nodes, node)
        node.clear_children()

    def write(self, outfile):
        """Write the graph in 'dot' format."""
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

    """Builder for Graph objects from gmxtree.GromacsTree representation."""

    def __init__(self, tree):
        """Initialize builder for a given tree representation."""
        self._tree = tree

    def _create_file_node(self, fileobj, filenodes):
        """Create graph node for a file object.

        filenodes is a dict() that maps file objects to their nodes, and is
        updated by this call.
        """
        nodename = re.subn(r'[-./]', '_', fileobj.get_relpath())[0]
        style = []
        properties = []
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
        node = Node(nodename, fileobj.get_name(), style, properties, is_file=True)
        filenodes[fileobj] = node
        return node

    def _get_file_edge_type(self, fromfile, tofile):
        """Get EdgeType for an edge between two file objects.

        Determines the type for the edge from the information provided by
        gmxtree.
        """
        intramodule = (fromfile.get_module() == tofile.get_module())
        is_legacy = not tofile.api_type_is_reliable()
        if fromfile.get_module() == tofile.get_module():
            return EdgeType.intramodule
        elif tofile.get_api_type() == DocType.internal and not tofile.is_public():
            if is_legacy:
                return EdgeType.legacy
            else:
                return EdgeType.undocumented
        elif fromfile.is_test_file():
            return EdgeType.test
        elif tofile.is_test_file():
            return EdgeType.undocumented
        elif fromfile.is_module_internal():
            if tofile.is_public():
                return EdgeType.pubimpl
            elif tofile.get_api_type() == DocType.library:
                return EdgeType.libimpl
            elif is_legacy:
                return EdgeType.legacy
            elif not tofile.is_documented():
                return EdgeType.undocumented
            else:
                raise ValueError('Unknown edge type between {0} and {1}'
                        .format(fromfile.get_relpath(), tofile.get_relpath()))
        elif fromfile.get_api_type() == DocType.library:
            return EdgeType.library
        elif fromfile.is_public() or fromfile.is_installed():
            if tofile.is_public() or tofile.is_installed():
                return EdgeType.public
            else:
                return EdgeType.undocumented
        elif is_legacy:
            return EdgeType.legacy
        else:
            raise ValueError('Unknown edge type between {0} and {1}'
                    .format(fromfile.get_relpath(), tofile.get_relpath()))

    def _create_file_edge(self, fromfile, tofile, filenodes):
        """Create edge between two file objects.

        Determines the type for the edge from the information provided by
        gmxtree.
        """
        edgetype = self._get_file_edge_type(fromfile, tofile)
        return Edge(filenodes[fromfile], filenodes[tofile], edgetype)

    def _create_file_edges(self, filenodes):
        """Create edges between all file nodes.

        Create edges between file nodes specified in filenodes from all include
        dependencies.  An edge is created only if both ends of the dependency
        are in the list of nodes.
        """
        edges = []
        for fileobj in filenodes.iterkeys():
            for includedfile in fileobj.get_includes():
                otherfile = includedfile.get_file()
                if otherfile and otherfile in filenodes:
                    edge = self._create_file_edge(fileobj, otherfile, filenodes)
                    edges.append(edge)
        return edges

    def _get_module_color(self, modulegroup):
        # If you change these styles, update also the legend in modulegraph.md
        if modulegroup == 'legacy':
            return 'fillcolor=grey75'
        elif modulegroup == 'analysismodules':
            return 'fillcolor="0 .2 1"'
        elif modulegroup == 'utilitymodules':
            return 'fillcolor=".08 .2 1"'
        elif modulegroup == 'mdrun':
            return 'fillcolor=".75 .2 1"'
        return None

    def _create_module_node(self, module):
        """Create node for a module."""
        style = []
        properties = []
        properties.append('shape=ellipse')
        if module.is_documented():
            properties.append('URL="\\ref {0}"'.format(module.get_name()))
        if not module.is_documented():
            fillcolor = self._get_module_color('legacy')
        else:
            fillcolor = self._get_module_color(module.get_group())
        if fillcolor:
            style.append('filled')
            properties.append(fillcolor)
        rootdir = module.get_root_dir()
        if rootdir.has_installed_files():
            properties.append('color=".66 .5 1"')
            properties.append('penwidth=3')
        nodename = 'module_' + re.subn(r'[-./]', '_', rootdir.get_relpath())[0]
        label = module.get_name()[7:]
        node = Node(nodename, label, style, properties)
        return node

    def _create_module_edges(self, modulenodes):
        """Create edges between all module nodes.

        Create edges between module nodes specified in modulenodes from all
        include dependencies.  An edge is created only if both ends of the
        dependency are in the list of nodes.
        """
        edges = []
        for moduleobj in modulenodes.iterkeys():
            for dep in moduleobj.get_dependencies():
                othermodule = dep.get_other_module()
                if othermodule and othermodule in modulenodes:
                    if dep.is_cycle_suppressed():
                        edgetype = EdgeType.cyclic
                    else:
                        edgetype = max([
                            self._get_file_edge_type(x.get_including_file(), x.get_file())
                            for x in dep.get_included_files()])
                    edge = Edge(modulenodes[moduleobj], modulenodes[othermodule], edgetype)
                    edges.append(edge)
        return edges

    def create_modules_graph(self):
        """Create module dependency graph."""
        nodes = []
        modulenodes = dict()
        libgromacsnode = Node('libgromacs', 'libgromacs')
        nodes.append(libgromacsnode)
        for moduleobj in self._tree.get_modules():
            node = self._create_module_node(moduleobj)
            if moduleobj.get_root_dir().get_relpath().startswith('src/gromacs'):
                libgromacsnode.add_child(node)
            else:
                nodes.append(node)
            modulenodes[moduleobj] = node
        edges = self._create_module_edges(modulenodes)
        graph = Graph(nodes, edges)
        graph.set_options(concentrate=False)
        return graph

    def create_module_file_graph(self, module):
        """Create file dependency graph for files within a module."""
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
    import os
    import sys

    from optparse import OptionParser

    from gmxtree import GromacsTree
    from reporter import Reporter

    parser = OptionParser()
    parser.add_option('-S', '--source-root',
                      help='Source tree root directory')
    parser.add_option('-B', '--build-root',
                      help='Build tree root directory')
    parser.add_option('--ignore-cycles',
                      help='Set file with module dependencies to ignore in cycles')
    parser.add_option('-o', '--outdir', default='.',
                      help='Specify output directory for graphs')
    parser.add_option('-q', '--quiet', action='store_true',
                      help='Do not write status messages')
    options, args = parser.parse_args()

    reporter = Reporter(quiet=True)

    if not options.quiet:
        sys.stderr.write('Scanning source tree...\n')
    tree = GromacsTree(options.source_root, options.build_root, reporter)
    tree.load_installed_file_list()
    if not options.quiet:
        sys.stderr.write('Reading source files...\n')
    tree.scan_files()
    if options.ignore_cycles:
        tree.load_cycle_suppression_list(options.ignore_cycles)
    if not options.quiet:
        sys.stderr.write('Reading Doxygen XML files...\n')
    tree.load_xml(only_files=True)

    if not options.quiet:
        sys.stderr.write('Writing graphs...\n')
    graphbuilder = GraphBuilder(tree)
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

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
