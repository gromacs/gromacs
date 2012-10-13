#!/usr/bin/python
"""Generate include dependency graphs for Gromacs.

This script can generate two types of include dependency graphs: per-file or
per-module (where module is equivalent to a subdirectory).
It is intended to be run on a subset of files under the src/ directory.
Output format is suitable for processing with 'dot'.

FILE GRAPHS

The legend for per-file graph nodex:
    gray:          source files
    light blue:    public headers
    dark blue:     library headers
    no background: other files

MODULE GRAPHS

Module graph will contain one node for each top-level subdirectory under src/,
except that the src/gromacs/ directory will be expanded one level further.

The legend for per-module graph links (a link with a certain color indicates
that types above it in the list are not present):
    red:          invalid dependency (e.g., undocumented file)
    dark blue:    library header depends on the other module
    light blue:   public header depends on the other module
    dashed black: source file depends on a library header in the other module
    solid black:  source file depends on a public header in the other module
    dotted grey:  test files depend on the other module
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

    priorities = {
            'unknown': 0,
            'undocumented': 1,
            'intramodule': 2,
            'library': 3,
            'public': 4,
            'libimpl': 5,
            'pubimpl': 6,
            'test': 7}

    def __init__(self, fromnode, tonode, link_type=None):
        self.fromnode = fromnode
        self.tonode = tonode
        self.link_type = link_type

    def refresh_type(self, reporter):
        """Initialize type of a link between two file nodes.

        Both endpoints of the link must be file objects when this method is
        called.
        """
        fromfile = self.fromnode.obj
        tofile = self.tonode.obj
        intramodule = \
                (fromfile.module.get_top_level_module() == \
                 tofile.module.get_top_level_module())
        if tofile.type != 'publicheader' and tofile.type != 'libheader':
            if intramodule:
                link_type = 'intramodule'
            else:
                reporter.error(fromfile.path,
                        'included file "{0}" is missing API definition'
                            .format(tofile.path))
                link_type = 'undocumented'
        elif fromfile.type == 'test':
            link_type = 'test'
        elif fromfile.type in ('source', 'header', 'implheader'):
            if tofile.type == 'publicheader':
                link_type = 'pubimpl'
            elif tofile.type == 'libheader':
                link_type = 'libimpl'
            else:
                reporter.error(fromfile.path,
                        'unknown link type to "{0}"'.format(tofile.path))
                link_type = 'unknown'
        elif fromfile.type == 'libheader':
            link_type = 'library'
        elif fromfile.type == 'publicheader':
            if tofile.type == 'publicheader' or tofile.doctype == 'public':
                link_type = 'public'
            else:
                reporter.error(fromfile.path,
                        'public API file includes non-public header "{0}"'
                            .format(tofile.path))
                link_type = 'undocumented'
        else:
            reporter.error(fromfile.path,
                    'unknown link type to "{0}"'.format(tofile.path))
            link_type = 'unknown'
        self.link_type = link_type

    def merge_link(self, other):
        """Merge another link into this one and choose an appropriate type.

        Updates the type of this link based on the types of the merged links.
        """
        if Link.priorities[other.link_type] < Link.priorities[self.link_type]:
            self.link_type = other.link_type

    def format(self):
        """Format this link for 'dot'."""
        if isinstance(self.fromnode.obj, File) and \
                isinstance(self.tonode.obj, File):
            properties = ''
        elif self.link_type == 'intramodule':
            properties = ''
        elif self.link_type == 'test':
            properties = 'color=grey75, style=dotted'
        elif self.link_type == 'libimpl':
            properties = 'color=".66 .5 1"'
        elif self.link_type == 'pubimpl':
            properties = 'color=".66 .2 1"'
        elif self.link_type == 'library':
            properties = 'color=black, style=dashed'
        elif self.link_type == 'public':
            properties = 'color=black'
        else: #unknown or undocumented
            properties = 'color=red'
        return '{0} -> {1} [{2}]'.format(self.fromnode.obj.nodename,
                                         self.tonode.obj.nodename,
                                         properties)

class Node(object):
    def __init__(self, obj):
        self.obj = obj
        self.children = []
        self.root = False

    def set_root(self):
        self.root = True

    def add_child(self, child):
        self.children.append(child)

    def clear_children(self):
        self.children = []

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
                              .format(self.obj.nodename)
                result += '        label = "{0}"\n'.format(self.obj.name)
            for child in self.children:
                result += child.format()
            if not self.root:
                result += '    }\n'
        else:
            result += '    {0} [{1}]\n'.format(
                    self.obj.nodename, self.obj.node_properties())
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


class File(object):
    def __init__(self, path, module):
        self.path = path
        self.name = os.path.basename(path)
        self.nodename = re.subn(r'[-./]', '_', path)[0]
        self.module = module
        if module.name == 'tests':
            self.type = 'test'
        elif re.search(r'\.c(pp)?$', self.name) != None:
            self.type = 'source'
        else:
            self.type = 'header'
        self.doctype = 'none'
        #headername = re.sub(r'\.cpp$', '.h', self.name)
        #implheadername = re.sub(r'\.cpp$', '-impl.h', self.name)
        self.links = []
        self.node = Node(self)
        self.installed = False

    def set_installed(self, reporter):
        if self.type != 'header':
            reporter.input_warning(self.path,
                    'installing {0} file'.format(self.type))
            return
        self.installed = True

    def add_dependency(self, dep):
        self.links.append(Link(self.node, dep.node))

    def get_node(self):
        return self.node

    def get_links(self):
        return self.links

    def node_properties(self):
        properties = []
        style = []
        properties.append('label="{0}"'.format(self.name))
        properties.append('URL="\\ref {0}"'.format(self.name))
        if not self.module:
            style.append('bold')
            properties.append('color=red')
        if self.type == 'source':
            style.append('filled')
            properties.append('fillcolor=grey75')
        elif self.type == 'publicheader':
            style.append('filled')
            properties.append('fillcolor=".66 .2 1"')
        elif self.type == 'libheader':
            style.append('filled')
            properties.append('fillcolor=".66 .5 1"')
        if style:
            properties.append('style="{0}"'.format(','.join(style)))
        return ', '.join(properties)

    def add_included_file(self, includefile, allfiles, includedirs, ignorelist,
            reporter):
        fullpath = os.path.join(includedirs[0], includefile)
        if os.path.abspath(fullpath) in ignorelist:
            return
        for includedir in includedirs:
            fullpath = os.path.abspath(os.path.join(includedir, includefile))
            if os.path.exists(fullpath):
                if self.installed and includedir != includedirs[0]:
                    reporter.error(self.path,
                            'installed header includes "{0}", '
                            'which is not found using relative path'
                                .format(includefile))
                includefile = fullpath
                break
        else:
            reporter.input_warning(self.path,
                    'included file "{0}" not found'
                        .format(includefile))
            return
        if includefile in allfiles:
            other = allfiles[includefile]
            if self.installed and not other.installed:
                reporter.error(self.path,
                        'installed header includes non-installed header "{0}"'
                            .format(other.path))
            self.add_dependency(other)
        #elif not dep in ignorelist:
        #    depfile = File(dep, None)
        #    files[dep] = depfile
        #    file.add_dependency(depfile)
        #    extrafiles.append(dep)

    def scan(self, filename, allfiles, includedirs, ignorelist, reporter):
        includedirs = [os.path.dirname(filename)] + includedirs
        infileblock = False
        foundfileblock = False
        docmodule = None
        with open(filename, 'r') as scanfile:
            for line in scanfile:
                match = re.match(r'#include "([^>"]*)"', line)
                if match:
                    self.add_included_file(match.group(1), allfiles,
                            includedirs, ignorelist, reporter)
                if not foundfileblock:
                    if infileblock:
                        if line.startswith(r' */'):
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
                            if docmodule:
                                reporter.error(self.path,
                                        'file documented in multiple modules')
                            docmodule = match.group(1)
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
        if self.doctype == 'none':
            reporter.error(self.path, 'file not documented')
        elif self.doctype == 'implementation' and \
                self.type in ('publicheader', 'libheader'):
            reporter.error(self.path,
                    'file documentation visibility incorrect')
        elif self.doctype == 'library' and self.type == 'publicheader':
            reporter.error(self.path,
                    'file documentation visibility incorrect')
        if self.installed and self.doctype not in ('public', 'unknown'):
            reporter.error(self.path,
                    'installed header has no public documentation')
        elif not self.installed and self.doctype == 'public':
            reporter.error(self.path,
                    'non-installed file has public documentation')
        selfmodnodename = self.module.nodename
        if docmodule and \
                not selfmodnodename.startswith('module_' + docmodule) and \
                not selfmodnodename.startswith('module_gromacs_' + docmodule):
            reporter.error(self.path,
                    'file documented in incorrect module "{0}"'
                        .format(docmodule))


class Module(object):
    def __init__(self, name, parent = None):
        self.parent = parent
        self.name = name
        if parent:
            self.nodename = parent.nodename + '_' + name
        else:
            self.nodename = 'module'
        self.files = []
        self.children = dict()
        self.is_top_level = (not parent or parent.name in ('', 'gromacs'))

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

    def create_node(self):
        node = Node(self)
        for childfile in self.files:
            node.add_child(childfile.get_node())
        for childmodule in self.children.itervalues():
            node.add_child(childmodule.create_node())
        return node

    def node_properties(self):
        properties = 'label="{0}", shape=ellipse'.format(self.name)
        properties += ', URL="\\ref module_{0}"'.format(self.name)
        return properties


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
        for scanfile in self.files.itervalues():
            for link in scanfile.get_links():
                link.refresh_type(reporter)

    def create_file_graph(self):
        rootnode = self.root.create_node()
        rootnode.set_root()
        links = []
        for scanfile in self.files.itervalues():
            links.extend(scanfile.get_links())
        graph = Graph([rootnode], links)
        return graph

    def create_modules_graph(self):
        rootnode = self.root.create_node()
        rootnode.set_root()
        links = []
        for scanfile in self.files.itervalues():
            links.extend(scanfile.get_links())
        graph = Graph([rootnode], links)
        for node in rootnode.get_children():
            if node.obj.name == 'gromacs':
                for child in node.get_children():
                    graph.collapse_node(child)
            else:
                graph.collapse_node(node)
        graph.set_options(concentrate=False)
        return graph

    def create_module_file_graph(self, module):
        rootnode = module.create_node()
        rootnode.set_root()
        links = []
        for scanfile in self.files.itervalues():
            links.extend(scanfile.get_links())
        graph = Graph([rootnode], links)
        graph.prune_links()
        return graph

    def get_toplevel_modules(self):
        result = []
        for module in self.root.children.itervalues():
            if module.name == 'gromacs':
                result.extend(module.children.itervalues())
            else:
                result.append(module)
        return result


def print_module_graph(outfile, deps, options):
    graph = deps.create_modules_graph()
    graph.write(outfile)

def print_file_graph(outfile, deps, options):
    graph = deps.create_file_graph()
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

def print_module_file_graph(outfile, deps, module, options):
    graph = deps.create_module_file_graph(module)
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
    parser.add_option('--left-to-right', action='store_true',
                      help='Lay out from left to right')
    parser.add_option('--file-graph',
                      help='Write graph for individual files')
    parser.add_option('--module-graph',
                      help='Write graph for modules')
    parser.add_option('--module-file-graphs', action='store_true',
                      help='Write file graphs for each module')
    options, args = parser.parse_args()

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
    if options.module_graph:
        graphpath = os.path.join(options.outdir, options.module_graph)
        with open(graphpath, 'w') as outfile:
            print_module_graph(outfile, deps, options)
    if options.file_graph:
        graphpath = os.path.join(options.outdir, options.file_graph)
        with open(graphpath, 'w') as outfile:
            print_file_graph(outfile, deps, options)
    if options.module_file_graphs:
        options.left_to_right = True
        for module in deps.get_toplevel_modules():
            filename = 'module_{0}-deps.dot'.format(module.name)
            with open(os.path.join(options.outdir, filename), 'w') as outfile:
                print_module_file_graph(outfile, deps, module, options)

main()
