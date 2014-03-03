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
import xml.etree.ElementTree as ET

def _show_list(title, objlist):
    if objlist:
        print '{0}:'.format(title)
        for obj in objlist:
            print '  ', obj

class Location(object):
    def __init__(self, elem):
        self.filepath = elem.attrib['file']
        self.line = int(elem.attrib['line'])
        self.column = elem.attrib['column']

    def __str__(self):
        return '{0}:{1}:{2}'.format(self.filepath, self.line, self.column)

class BodyLocation(object):
    def __init__(self, elem):
        self.filepath = elem.attrib['bodyfile']
        self.startline = int(elem.attrib['bodystart'])
        self.endline = int(elem.attrib['bodyend'])

    def __cmp__(self, other):
        result = cmp(self.filepath, other.filepath)
        if result == 0:
            result = cmp(self.startline, other.startline)
        if result == 0:
            result = cmp(self.endline, other.endline)
        return result

    def __hash__(self):
        return hash(self.filepath) ^ hash(self.startline) ^ hash(self.endline)

    def __str__(self):
        return '{0}:{1}-{2}'.format(self.filepath, self.startline, self.endline)

class LocationWithBody(object):
    def __init__(self, elem):
        self._location = Location(elem)
        if 'bodyfile' in elem.attrib:
            self._bodylocation = BodyLocation(elem)
        else:
            self._bodylocation = None

    def __str__(self):
        if not self._bodylocation:
            return '{0} (no body)'.format(self._location)
        else:
            return '{0} / {1}'.format(self._location, self._bodylocation)

    def get_location(self):
        return self._location

    def get_body_location(self):
        return self._bodylocation

    def has_same_body_location(self):
        return self._location.filepath == self._bodylocation.filepath and \
                self._location.line == self._bodylocation.startline

class MemberSection(object):
    def __init__(self, kind):
        self._kind = kind
        self._members = []

    def __str__(self):
        return self._kind

    def add_member(self, member):
        self._members.append(member)

    def replace_member(self, old, new):
        try:
            pos = self._members.index(old)
        except ValueError:
            return
        self._members[pos] = new

class Entity(object):
    def __init__(self, name, refid):
        self._name = name
        self._id = refid
        self._has_brief_description = False
        self._has_detailed_description = False
        self._has_inbody_description = False
        self._visibility = 'none'

    def __str__(self):
        return self._name

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_visibility(self):
        return self._visibility

    def is_documented(self):
        return self._visibility != 'none'

    def has_brief_description(self):
        return self._has_brief_description

    def has_inbody_description(self):
        return self._has_inbody_description

    def _process_descriptions(self, briefelem, detailselem, inbodyelem):
        if briefelem is not None and len(briefelem) > 0:
            self._has_brief_description = True
            self._visibility = 'public'
        if detailselem is not None and len(detailselem) > 0:
            self._visibility = 'public'
            # Gromacs-specific:
            # \internal is used at the beginning of a comment block to
            # mark the block internal to the module.
            # \libinternal is used similarly, and inserts custom XML
            # elements.
            if detailselem[0].tag == 'internal':
                if len(detailselem) > 1:
                    sys.stderr.write('warning: unexpected extra elements after <internal>\n')
                self._visibility = 'internal'
            if detailselem[0].find('libinternal') is not None:
                if self._visibility != 'public':
                    sys.stderr.write('warning: multiple visibility specifications\n')
                self._visibility = 'library'
            self._has_detailed_description = True
        if inbodyelem is not None:
            self._has_inbody_description = (len(inbodyelem) > 0)

    def show_base(self):
        print 'ID:         {0}'.format(self._id)
        print 'Name:       {0}'.format(self._name)
        doctype = []
        if self._has_brief_description:
            doctype.append('brief')
        if self._has_detailed_description:
            doctype.append('details')
        if self._has_inbody_description:
            doctype.append('in-body')
        if not doctype:
            doctype.append('none')
        print 'Doc:        {0}'.format(', '.join(doctype))
        print 'Visibility: {0}'.format(self._visibility)

class Compound(Entity):
    def __init__(self, name, refid):
        Entity.__init__(self, name, refid)
        self._docset = None
        self._members = dict()
        self._children = set()
        self._sections = []
        self._groups = set()
        self._loaded = False

    def set_documentation_set(self, docset):
        self._docset = docset

    def add_member(self, member):
        self._members[member.get_id()] = member

    def add_group(self, compound):
        self._groups.add(compound)

    def replace_member(self, old, new):
        if old.get_id() not in self._members:
            raise ValueError("Trying to replace a non-existent member")
        elif new.get_id() in self._members:
            raise ValueError("Trying to replace with an existing member")
        self._members[old.get_id()] = new
        for section in self._sections:
            section.replace_member(old, new)

    def get_xml_path(self):
        return os.path.join(self._docset.get_xmlroot(), self.get_id() + '.xml')

    def get_groups(self):
        return self._groups

    def load_details(self):
        if self._loaded:
            return
        xmlfile = self.get_xml_path()
        compoundtree = ET.parse(xmlfile)
        root = compoundtree.getroot()
        if len(root) > 1:
            sys.stderr.write("warning: more than one compound in a file\n")
        briefelem = None
        detailselem = None
        missing_members = set(self._members.values())
        for elem in root.find('compounddef'):
            if elem.tag == 'compoundname':
                if elem.text != self.get_name():
                    sys.stderr.write("warning: compound name mismatch: '{0}' vs '{1}'\n".format(self.get_name(), elem.text))
            elif elem.tag == 'briefdescription':
                briefelem = elem
            elif elem.tag == 'detaileddescription':
                detailselem = elem
            elif elem.tag in ('includes', 'includedby', 'incdepgraph',
                    'invincdepgraph', 'inheritancegraph', 'collaborationgraph',
                    'programlisting', 'templateparamlist', 'listofallmembers'):
                pass
            elif elem.tag.startswith('inner'):
                refid = elem.attrib['refid']
                reftype = elem.tag[5:]
                # TODO: Handle 'prot' attribute?
                refcompound = self._docset.get_compound(refid)
                self._children.add(refcompound)
                if reftype == 'file':
                    self._load_inner_file(refcompound)
                elif reftype == 'dir':
                    self._load_inner_dir(refcompound)
                elif reftype == 'group':
                    self._load_inner_group(refcompound)
                elif reftype == 'namespace':
                    self._load_inner_namespace(refcompound)
                elif reftype == 'class':
                    self._load_inner_class(refcompound)
                else:
                    sys.stderr.write("warning: unknown inner compound type '{0}'\n".format(reftype))
            elif elem.tag == 'sectiondef':
                # TODO: Handle header and description elements
                kind = elem.attrib['kind']
                section = MemberSection(kind)
                self._sections.append(section)
                for memberelem in elem.iter('memberdef'):
                    refid = memberelem.attrib['id']
                    member = self._members[refid]
                    member.load_details_from_element(memberelem)
                    section.add_member(member)
                    if member in missing_members:
                        missing_members.remove(member)
            else:
                if not self._load_element(elem):
                    sys.stderr.write("warning: unknown compound child element '{0}'\n".format(elem.tag))
        # TODO: Test again once enum values are interpreted
        #if missing_members:
        #    sys.stderr.write('warning: members without section\n')
        self._process_descriptions(briefelem, detailselem, None)
        self._loaded = True

    def _load_inner_file(self, compound):
        sys.stderr.write("warning: unexpected inner file\n")

    def _load_inner_dir(self, compound):
        sys.stderr.write("warning: unexpected inner dir\n")

    def _load_inner_group(self, compound):
        sys.stderr.write("warning: unexpected inner group\n")

    def _load_inner_namespace(self, compound):
        sys.stderr.write("warning: unexpected inner namespace\n")

    def _load_inner_class(self, compound):
        sys.stderr.write("warning: unexpected inner class\n")

    def _load_element(self, element):
        return False

    def _ensure_loaded(self):
        self.load_details()

    def show_base(self):
        Entity.show_base(self)
        if self._groups:
            print 'Groups:   {0}'.format(', '.join(map(str, self._groups)))

    def show_members(self):
        for section in self._sections:
            print 'Member section: {0}'.format(section)
            for member in section._members:
                print '  ', member

class File(Compound):
    def __init__(self, name, refid):
        Compound.__init__(self, name, refid)
        self._path = None
        self._directory = None
        self._classes = set()
        self._namespaces = set()

    def _load_inner_class(self, compound):
        compound.add_file(self)
        self._classes.add(compound)

    def _load_inner_namespace(self, compound):
        compound.add_file(self)
        self._namespaces.add(compound)

    def _load_element(self, elem):
        if elem.tag == 'location':
            self._path = elem.attrib['file']
            return True
        return False

    def set_directory(self, directory):
        self._directory = directory

    def get_path(self):
        return self._path

    def get_directory(self):
        return self._directory

    def show(self):
        self._ensure_loaded()
        self.show_base()
        print 'Path:      {0}'.format(self._path)
        print 'Directory: {0}'.format(self._directory)
        _show_list('Namespaces', self._namespaces)
        _show_list('Classes', self._classes)
        self.show_members()

class Directory(Compound):
    def __init__(self, name, refid):
        Compound.__init__(self, name, refid)
        self._path = None
        self._parent = None
        self._subdirs = set()
        self._files = set()

    def _load_inner_file(self, compound):
        compound.set_directory(self)
        self._files.add(compound)

    def _load_inner_dir(self, compound):
        compound._parent = self
        self._subdirs.add(compound)

    def _load_element(self, elem):
        if elem.tag == 'location':
            self._path = elem.attrib['file']
            return True
        return False

    def get_path(self):
        return self._path

    def get_parent(self):
        return self._parent

    def get_subdirectories(self):
        return self._subdirs

    def show(self):
        self._ensure_loaded()
        self.show_base()
        print 'Path:      {0}'.format(self._path)
        if self._parent:
            print 'Parent:    {0}'.format(self._parent)
        _show_list('Subdirectories', self._subdirs)
        _show_list('Files', self._files)

class Group(Compound):
    def __init__(self, name, refid):
        Compound.__init__(self, name, refid)
        self._title = None
        self._files = set()
        self._nestedgroups = set()
        self._namespaces = set()
        self._classes = set()

    def _load_inner_file(self, compound):
        compound.add_group(self)
        self._files.add(compound)

    # Doxygen 1.8.5 doesn't seem to put the directories into the XML output,
    # even though they are in the HTML output as group members...

    def _load_inner_group(self, compound):
        compound.add_group(self)
        self._nestedgroups.add(compound)

    def _load_inner_namespace(self, compound):
        compound.add_group(self)
        self._namespaces.add(compound)

    def _load_inner_class(self, compound):
        compound.add_group(self)
        self._classes.add(compound)

    def _load_element(self, elem):
        if elem.tag == 'title':
            self._title = elem.text
            return True
        return False

    def show(self):
        self._ensure_loaded()
        self.show_base()
        print 'Title:     {0}'.format(self._title)
        print 'Inner compounds:'
        for compound in self._children:
            print '  ', compound
        self.show_members()

class Namespace(Compound):
    def __init__(self, name, refid):
        Compound.__init__(self, name, refid)
        self._doclocation = None
        self._files = set()
        self._parent = None
        self._innernamespaces = set()
        self._classes = set()

    def _load_inner_namespace(self, compound):
        compound._parent = self
        self._innernamespaces.add(compound)

    def _load_inner_class(self, compound):
        compound.set_namespace(self)
        self._classes.add(compound)

    def _load_element(self, elem):
        if elem.tag == 'location':
            self._doclocation = Location(elem)
            return True
        return False

    def add_file(self, compound):
        self._files.add(compound)

    def get_location(self):
        return self._doclocation

    def show(self):
        self._ensure_loaded()
        self.show_base()
        print 'Doc. loc.: {0}'.format(self._doclocation)
        _show_list('Inner namespaces', self._innernamespaces)
        _show_list('Classes', self._classes)
        self.show_members()

class Class(Compound):
    def __init__(self, name, refid):
        Compound.__init__(self, name, refid)
        self._location = None
        self._namespace = None
        self._files = set()
        self._baseclasses = []
        self._derivedclasses = set()
        self._outerclass = None
        self._innerclasses = set()

    def _load_inner_class(self, compound):
        compound.set_outer_class(self)
        self._innerclasses.add(compound)

    def _load_element(self, elem):
        if elem.tag == 'basecompoundref':
            # TODO: Handle unknown bases?
            if 'refid' in elem.attrib:
                refid = elem.attrib['refid']
                # TODO: Handle prot and virt attributes, check name?
                base = self._docset.get_compound(refid)
                self._baseclasses.append(base)
            return True
        if elem.tag == 'derivedcompoundref':
            refid = elem.attrib['refid']
            # TODO: Handle prot and virt attributes, check name?
            derived = self._docset.get_compound(refid)
            self._derivedclasses.add(derived)
            return True
        elif elem.tag == 'location':
            self._location = LocationWithBody(elem)
            return True
        return False

    def add_file(self, compound):
        self._files.add(compound)

    def set_namespace(self, compound):
        self._namespace = compound

    def set_outer_class(self, compound):
        self._outerclass = compound

    def get_location(self):
        return self._location.get_location()

    def show(self):
        self._ensure_loaded()
        self.show_base()
        print 'Location:   {0}'.format(self._location)
        print 'Namespace:  {0}'.format(self._namespace)
        if self._outerclass:
            print 'Outer cls:  {0}'.format(self._outerclass)
        _show_list('Inner classes', self._innerclasses)
        self.show_members()

class Member(Entity):
    def __init__(self, name, refid):
        Entity.__init__(self, name, refid)
        self._parents = set()
        self._location = None
        self._alternates = set()
        self._loaded = False

    def add_parent_compound(self, compound):
        self._parents.add(compound)

    def get_parent_compounds(self):
        return self._parents

    def get_inherited_visibility(self):
        result = 'none'
        for parent in self._parents:
            if parent.get_visibility() == 'public':
                result = 'public'
                break
            elif parent.get_visibility() == 'library':
                if result in ('none', 'internal'):
                    result = 'library'
            elif parent.get_visibility() == 'internal':
                if result == 'none':
                    result = 'internal'
        return result

    def has_same_body_location(self):
        return self._location.has_same_body_location()

    def get_location(self):
        return self._location.get_location()

    def get_body_location(self):
        # TODO: Check whether this can be removed once enum values are handled
        if not self._location:
            return None
        return self._location.get_body_location()

    def merge_definition(self, definition):
        self._parents.update(definition._parents)
        self._alternates.add(definition)

    def load_details_from_element(self, rootelem):
        if self._loaded:
            # TODO: It would be nice to verify that the same information
            # is present in all instances
            return
        # TODO: Process the attributes
        briefelem = None
        detailselem = None
        inbodyelem = None
        for elem in rootelem:
            if elem.tag == 'name':
                if elem.text != self.get_name():
                    sys.stderr.write("warning: member name mismatch: '{0}' vs '{1}'\n".format(self.get_name(), elem.text))
            elif elem.tag == 'briefdescription':
                briefelem = elem
            elif elem.tag == 'detaileddescription':
                detailselem = elem
            elif elem.tag == 'inbodydescription':
                inbodyelem = elem
            elif elem.tag == 'location':
                self._location = LocationWithBody(elem)
            else:
                # TODO Process the rest of the elements
                pass
        self._process_descriptions(briefelem, detailselem, inbodyelem)
        self._loaded = True

    def show(self):
        self.show_base()
        print 'Parent vis: {0}'.format(self.get_inherited_visibility())
        print 'Location:   {0}'.format(self._location)
        _show_list('Parents', self._parents)

class DocumentationSet(object):
    def __init__(self, xmlroot):
        self._xmlroot = xmlroot
        indextree = ET.parse(os.path.join(xmlroot, 'index.xml'))
        self._compounds = dict()
        self._members = dict()
        self._files = dict()
        for compoundelem in indextree.getroot():
            name = compoundelem.find('name').text
            refid = compoundelem.attrib['refid']
            kind = compoundelem.attrib['kind']
            if kind == 'file':
                compound = File(name, refid)
            elif kind == 'dir':
                compound = Directory(name, refid)
            elif kind == 'group':
                compound = Group(name, refid)
            elif kind == 'namespace':
                compound = Namespace(name, refid)
            elif kind in ('class', 'struct', 'union'):
                compound = Class(name, refid)
            elif kind in ('page', 'example'):
                continue
            else:
                sys.stderr.write("warning: unknown compound kind '{0}'\n".format(kind))
                continue
            compound.set_documentation_set(self)
            self._compounds[refid] = compound
            for memberelem in compoundelem.iter('member'):
                name = memberelem.find('name').text
                refid = memberelem.attrib['refid']
                kind = memberelem.attrib['kind']
                if refid in self._members:
                    member = self._members[refid]
                else:
                    member = Member(name, refid)
                    self._members[refid] = member
                member.add_parent_compound(compound)
                compound.add_member(member)

    def load_details(self):
        for compound in self._compounds.itervalues():
            compound.load_details()
            if isinstance(compound, File):
                self._files[compound.get_path()] = compound
        # TODO: Add links to files using location

    def remove_duplicates(self):
        members_by_body = dict()
        for member in self._members.itervalues():
            bodyloc = member.get_body_location()
            if bodyloc:
                index = (bodyloc, member.get_name())
                if index not in members_by_body:
                    members_by_body[index] = []
                members_by_body[index].append(member)
        for memberlist in members_by_body.itervalues():
            if len(memberlist) > 1:
                if len(memberlist) > 2:
                    # TODO: Remove duplicate prototypes to avoid these warnings
                    sys.stderr.write("warning: more than two duplicates for a member '{0}'\n".format(memberlist[0]))
                    continue
                if memberlist[0].has_same_body_location():
                    declaration = memberlist[1]
                    definition = memberlist[0]
                elif memberlist[1].has_same_body_location():
                    declaration = memberlist[0]
                    definition = memberlist[1]
                else:
                    # TODO: Remove duplicate prototypes to avoid these warnings
                    sys.stderr.write("warning: duplicate member definition '{0}', but neither matches the body location\n".format(memberlist[0]))
                    continue
                self._members[definition.get_id()] = declaration
                declaration.merge_definition(definition)
                for compound in definition.get_parent_compounds():
                    compound.replace_member(definition, declaration)

    def get_xmlroot(self):
        return self._xmlroot

    def get_compound(self, refid):
        return self._compounds[refid]

    def get_member(self, refid):
        return self._members[refid]

    def get_compounds(self, types, predicate=None):
        result = []
        for compound in self._compounds.itervalues():
            if isinstance(compound, types) and \
                    (predicate is None or predicate(compound)):
                result.append(compound)
        return result

    def get_members(self, types=None, predicate=None):
        # self._members can contain duplicates
        result = set()
        for member in self._members.itervalues():
            if (types is None or isinstance(member, types)) and \
                    (predicate is None or predicate(member)):
                result.add(member)
        return list(result)

    def get_files(self, paths=None):
        if paths:
            return self.get_compounds(File, lambda x: x.get_name().endswith(paths))
        else:
            return self.get_compounds(File)

    def get_directories(self, paths):
        return self.get_compounds(Directory, lambda x: x.get_name().endswith(paths))

    def get_groups(self, name):
        return self.get_compounds(Group, lambda x: x.get_name() in name)

    def get_namespaces(self, name):
        return self.get_compounds(Namespace, lambda x: x.get_name() in name)

    def get_classes(self, name):
        return self.get_compounds(Class, lambda x: x.get_name() in name)

    def get_functions(self, name):
        return self.get_members(Member, lambda x: x.get_name() in name)

def main():
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-R', '--root-dir',
                      help='Doxygen XML root directory')
    parser.add_option('-F', '--show-file', action='append',
                      help='Show contents of given file')
    parser.add_option('-d', '--show-dir', action='append',
                      help='Show contents of given directory')
    parser.add_option('-g', '--show-group', action='append',
                      help='Show contents of given group')
    parser.add_option('-n', '--show-namespace', action='append',
                      help='Show contents of given namespace')
    parser.add_option('-c', '--show-class', action='append',
                      help='Show contents of given class')
    parser.add_option('-f', '--show-function', action='append',
                      help='Show details of given function')
    options, args = parser.parse_args()

    sys.stderr.write('Loading index.xml...\n')
    docset = DocumentationSet(options.root_dir)
    sys.stderr.write('Loading details...\n')
    docset.load_details()
    sys.stderr.write('Processing...\n')
    docset.remove_duplicates()
    objlist = []
    if options.show_file:
        objlist.extend(docset.get_files(tuple(options.show_file)))
    if options.show_dir:
        objlist.extend(docset.get_directories(tuple(options.show_dir)))
    if options.show_group:
        objlist.extend(docset.get_groups(tuple(options.show_group)))
    if options.show_namespace:
        # TODO: Replace file names with anonymous_namespace{filename}
        objlist.extend(docset.get_namespaces(tuple(options.show_namespace)))
    if options.show_class:
        objlist.extend(docset.get_classes(tuple(options.show_class)))
    if options.show_function:
        objlist.extend(docset.get_functions(tuple(options.show_function)))
    for obj in objlist:
        obj.show()
if __name__ == '__main__':
    main()
