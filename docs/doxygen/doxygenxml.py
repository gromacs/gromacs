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

"""Doxygen XML output parser.

This module implements a parser for the Doxygen XML output, converting it into
an object model that can be used to navigate the documentation.  It also uses
knowledge from how Doxygen works to provide access to things like visibility of
individual member documentation (e.g., based on what is the visibility of its
parent compound objects).

The object model is rooted at a DocumentationSet object.  Each documented
entity is modeled as an Entity, and this has subclasses Member and Compound to
correspond to the two categories of items that Doxygen handles.  These classes
are further subclassed to match each kind of entity that Doxygen produces.
Only kinds produced by Doxygen from C/C++ code are modeled.  Everything else
is ignored after a warning.

Currently the member entities are not completely parsed from the XML files, and
the interface may need additional work to provide convenient access to all
member types and their common properties.  For now, focus is in modeling the
compound entities.

The implementation is mostly independent of any GROMACS-specific rules, except
for the following:
 - DocType.library is a GROMACS-specific construct that is deduced from the
   contents of the detailed description (presence of a \libinternal command in
   the Doxygen comment triggers it).
 - DocType.internal is deduced from the presence of a \internal command that
   covers the whole detailed description.
 - List of extensions for determining whether a file is a source file only
   contains extensions actually used by GROMACS.
It would be possible to move these out from this file, but that would require
exposing the XML representation for the descriptions, which is not nice either.

The module can also be run as a script that can dump out different parts of the
object model.  This can be used to debug the parser, as well as check what is
actually in the XML documentation.
"""

import os.path
import xml.etree.ElementTree as ET

import reporter

#####################################################################
# Helper functions and classes

def _show_list(title, objlist):
    """Helper function for formatting a list of objects for debug output."""
    if objlist:
        print '{0}:'.format(title)
        for obj in objlist:
            print '  ', obj

class DocType(object):

    """Documentation visibility in the generated documentation."""

    # Mapping to string representations for the internal integer values
    _names = ['undocumented', 'internal', 'library', 'public']

    def __init__(self, value):
        """Initialize a DocType instance.

        DocType.{none,internal,library,public} should be used outside the class
        instead of calling the constructor.
        """
        self._value = value

    def __str__(self):
        """Return string representation for the documentation type."""
        return self._names[self._value]

    def __cmp__(self, other):
        """Order documentation types in the order of visibility."""
        return cmp(self._value, other._value)

# Static values for documentation types.
DocType.none = DocType(0)
DocType.internal = DocType(1)
DocType.library = DocType(2)
DocType.public = DocType(3)

class Location(object):

    """Location of a Doxygen entity.

    This class contains the logic to parse a <location> tag in Doxygen XML.
    It is used as the entity location in cases where body location is not
    expected, or as part of a LocationWithBody.
    """

    def __init__(self, elem):
        """Initialize location from a <location> element."""
        self.filepath = elem.attrib['file']
        self.line = int(elem.attrib['line'])
        self.column = elem.attrib['column']

    def __str__(self):
        return '{0}:{1}'.format(self.filepath, self.line)

    def get_reporter_location(self):
        return reporter.Location(self.filepath, self.line)

    def get_full_string(self):
        return '{0}:{1}:{2}'.format(self.filepath, self.line, self.column)

class BodyLocation(object):

    """Body location of a Doxygen entity.

    This class contains the logic to parse a body location from a <location>
    tag in Doxygen XML.  Not all entities have these attributes.
    This is only used as part of a LocationWithBody, which handles cases where
    the body location is optional.

    The body location can be compared and hashed so that it can be used in
    a dictionary for DocumentationSet.merge_duplicates().
    """

    def __init__(self, elem):
        """Initialize body location from a <location> element."""
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
        return '{0}:{1}'.format(self.filepath, self.startline)

    def get_full_string(self):
        if self.endline < 0:
            return self.__str__()
        return '{0}:{1}-{2}'.format(self.filepath, self.startline, self.endline)

class LocationWithBody(object):

    """Location for a Doxygen entity that can have a body location.

    This class is used to represent the location of a Doxygen entity that can
    have a body location.
    """

    def __init__(self, elem):
        """Initialize location from a <location> element."""
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

    def get_reporter_location(self):
        """Return reporter location for this location.

        All issues are reported at the main location, which should match with
        the declaration, where most of the documentation typically is.
        """
        return self._location.get_reporter_location()

    def get_location(self):
        return self._location

    def get_body_location(self):
        return self._bodylocation

    def has_same_body_location(self):
        """Check whether main location matches body location.

        If the main location is different, then it likely points to the
        declaration of the function.
        """
        return self._location.filepath == self._bodylocation.filepath and \
                self._location.line == self._bodylocation.startline

class MemberSection(object):

    """Section of members within a compound entity."""

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

#####################################################################
# Documentation entities

class Entity(object):

    """Doxygen documentation entity.

    This class represents common properties of an entity that can contain
    Doxygen documentation.
    """

    def __init__(self, name, refid):
        self._docset = None
        self._name = name
        self._id = refid
        self._has_brief_description = False
        self._has_detailed_description = False
        self._has_inbody_description = False
        self._visibility = DocType.none

    def __str__(self):
        return self._name

    def _get_reporter(self):
        """Return reporter to use for parsing issues."""
        return self._docset.get_reporter()

    def set_documentation_set(self, docset):
        """Set the documentation set this entity belongs to.

        The documentation set parent provides access to a common reporter
        object, and also allows the entity to resolve references to other
        entities while loading XML information.
        """
        assert self._docset is None
        self._docset = docset

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_reporter_location(self):
        return reporter.Location('<{0}>'.format(self._name), None)

    def get_visibility(self):
        return self._visibility

    def is_documented(self):
        return self._visibility != DocType.none

    def has_brief_description(self):
        return self._has_brief_description

    def has_inbody_description(self):
        return self._has_inbody_description

    def _process_descriptions(self, briefelem, detailselem, inbodyelem):
        reporter = self._get_reporter()
        if briefelem is not None and len(briefelem) > 0:
            self._has_brief_description = True
            self._visibility = DocType.public
        if detailselem is not None and len(detailselem) > 0:
            self._visibility = DocType.public
            # Gromacs-specific:
            # \internal is used at the beginning of a comment block to
            # mark the block internal to the module.
            # \libinternal is used similarly, and inserts custom XML
            # elements.
            if detailselem[0].tag == 'internal':
                if len(detailselem) == 1:
                    self._visibility = DocType.internal
                else:
                    # TODO: Should we also check if internal appears elsewhere?
                    reporter.doc_note(self, '\internal does not cover whole documentation')
            if detailselem[0].find('libinternal') is not None:
                if self._visibility == DocType.public:
                    self._visibility = DocType.library
                else:
                    reporter.doc_error(self, '\libinternal should not be used inside \internal')
            self._has_detailed_description = True
        if inbodyelem is not None:
            self._has_inbody_description = (len(inbodyelem) > 0)

    def show_base(self):
        """Format information for common properties.

        This is called from subclass show() methods to show base information
        about the entity.
        """
        print 'ID:         {0}'.format(self._id)
        print 'Name:       {0}'.format(self._name)
        print 'Location:   {0}'.format(self.get_reporter_location())
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

# Member entities

class Member(Entity):

    """Member entity.

    In Doxygen, a member entity is an entity such as a function or an enum that
    cannot contain other documented entities (an enum is a slight exception, as
    enum values are still nested within the enum member).  A member always
    belongs to one (or more) compounds, which means that the detailed
    documentation for the member appears on the documentation page for that
    compound.  If none of the parent compounds are documented, the member
    doesn't appear anywhere, even if it is documented.

    Member information is loaded from a parent compound's XML file.  If there
    is more than one parent, the first one encountered will be used
    (presumably, Doxygen duplicates the information into each XML file).
    """

    def __init__(self, name, refid):
        Entity.__init__(self, name, refid)
        self._parents = set()
        self._class = None
        self._namespace = None
        self._files = set()
        self._group = None
        self._location = None
        self._alternates = set()
        self._loaded = False
        # TODO: Move to Entity?
        self._xmlpath = None

    def add_parent_compound(self, compound):
        """Add a compound that contains this member."""
        self._parents.add(compound)
        if isinstance(compound, Class):
            assert self._class is None
            self._class = compound
        elif isinstance(compound, Namespace):
            assert self._namespace is None
            self._namespace = compound
        elif isinstance(compound, File):
            self._files.add(compound)
        elif isinstance(compound, Group):
            assert self._group is None
            self._group = compound
        else:
            assert False

    def merge_definition(self, definition):
        """Merge another member into this.

        See DocumentationSet.merge_duplicates().
        """
        assert self._class is None
        assert definition._class is None
        assert self._group == definition._group
        assert self._namespace == definition._namespace
        self._parents.update(definition._parents)
        self._files.update(definition._files)
        self._alternates.add(definition)

    def load_details_from_element(self, rootelem, xmlpath):
        """Load details for the member from a given XML element.

        This method is called when encountering member definitions while
        processing a compound XML file to load the information for that member.
        It processes common properties for a member, and delegates other
        elements to _load_element().
        """
        if self._loaded:
            # TODO: It would be nice to verify that the same information
            # is present in all instances
            return
        self._xmlpath = xmlpath
        # TODO: Process the attributes
        reporter = self._get_reporter()
        briefelem = None
        detailselem = None
        inbodyelem = None
        for elem in rootelem:
            if elem.tag == 'name':
                if elem.text != self.get_name():
                    reporter.xml_assert(xmlpath,
                            "member name mismatch: '{0}' (in index.xml) vs. '{1}'".format(
                                self.get_name(), elem.text))
            elif elem.tag == 'briefdescription':
                briefelem = elem
            elif elem.tag == 'detaileddescription':
                detailselem = elem
            elif elem.tag == 'inbodydescription':
                # TODO: in-body description is probably only possible for
                # functions; move it there.
                inbodyelem = elem
            elif elem.tag == 'location':
                self._location = LocationWithBody(elem)
            else:
                if not self._load_element(elem):
                    # TODO Process the rest of the elements so that we can check this
                    #reporter.xml_assert(xmlpath,
                    #        "unknown member child element '{0}'".format(elem.tag))
                    pass
        self._process_descriptions(briefelem, detailselem, inbodyelem)
        self._loaded = True

    def _load_element(self, element):
        """Load data from a child XML element.

        This method is called for all XML elements under the <memberdef>
        element that are not handled directly by the Member class.
        Derived classes should return True if they process the element.
        """
        return False

    def _get_raw_location(self):
        """Returns the BodyLocation object associated with this member.

        This is necessary so that EnumValue can override it report a non-empty
        location: Doxygen doesn't provide any location for <enumvalue>.
        """
        return self._location

    def get_reporter_location(self):
        return self._get_raw_location().get_reporter_location()

    def get_location(self):
        """Return main location for the member.

        This typically corresponds to the declaration.
        """
        return self._get_raw_location().get_location()

    def get_body_location(self):
        """Return location of the body for the member.

        Some types of members do not have a body location, in which case this
        returns None.
        """
        return self._get_raw_location().get_body_location()

    def has_same_body_location(self):
        """Check whether the main location is the same as body location."""
        return self._get_raw_location().has_same_body_location()

    def get_namespace(self):
        return self._namespace

    def get_parent_compounds(self):
        return self._parents

    def get_inherited_visibility(self):
        return max([parent.get_visibility() for parent in self._parents])

    def show(self):
        self.show_base()
        if self._alternates:
            idlist = [x.get_id() for x in self._alternates]
            print 'Alt. IDs:   {0}'.format(', '.join(idlist))
        print 'Parent vis: {0}'.format(self.get_inherited_visibility())
        print 'Location:   {0}'.format(self.get_location().get_full_string())
        print 'Body loc:   {0}'.format(self.get_body_location().get_full_string())
        _show_list('Parents', self._parents)

class Define(Member):
    pass

class Variable(Member):
    pass

class Typedef(Member):
    pass

class Enum(Member):
    def __init__(self, name, refid):
        Member.__init__(self, name, refid)
        self._values = set()

    def _load_element(self, elem):
        if elem.tag == 'enumvalue':
            refid = elem.attrib['id']
            # Doxygen seems to sometimes assign the same ID to a singleton enum
            # value (this already triggers a warning in loading index.xml).
            if refid == self.get_id():
                return True
            member = self._docset.get_member(refid)
            member.set_enum(self)
            member.load_details_from_element(elem, self._xmlpath)
            self._values.add(member)
            return True
        return False

    def get_values(self):
        return self._values

class EnumValue(Member):
    def __init__(self, name, refid):
        Member.__init__(self, name, refid)
        self._enum = None

    def set_enum(self, member):
        assert self._enum is None
        self._enum = member

    def _get_raw_location(self):
        return self._enum._get_raw_location()

class Function(Member):
    pass

class FriendDeclaration(Member):
    pass

# Compound entities

class Compound(Entity):

    """Compound entity.

    In Doxygen, a compound entity is an entity that has its own documentation
    page, and can contain other documented entities (either members, or other
    compounds).  Examples of compounds are files and classes.
    A compound entity always appears in the documentation, even if it is
    contained in another compound that is not documented.

    The list of members for a compound is initialized when the XML index file
    is read.  All other information is loaded from an XML file that is specific
    to the compound.  In addition to describing the compound, this XML file
    contains references to contained compounds, and details of all members
    within the compound.
    """

    def __init__(self, name, refid):
        Entity.__init__(self, name, refid)
        self._members = dict()
        self._children = set()
        self._sections = []
        self._groups = set()
        self._loaded = False

    def _get_xml_path(self):
        """Return path to the details XML file for this compound."""
        return os.path.join(self._docset.get_xmlroot(), self.get_id() + '.xml')

    def add_member(self, member):
        """Add a contained member."""
        self._members[member.get_id()] = member

    def add_group(self, compound):
        """Add a group (a compound entity) that contains this entity."""
        self._groups.add(compound)

    def replace_member(self, old, new):
        if old.get_id() not in self._members:
            raise ValueError("Trying to replace a non-existent member")
        elif new.get_id() in self._members:
            raise ValueError("Trying to replace with an existing member")
        self._members[old.get_id()] = new
        for section in self._sections:
            section.replace_member(old, new)

    def load_details(self):
        """Load details for the compound from its details XML file.

        This method processes common properties for a compound.
        References to inner compounds are delegated to _load_inner_*() methods,
        and all members encountered in the XML file are loaded with
        Member.load_details_from_element().
        Other elements are delegated to _load_element().
        """
        if self._loaded:
            return
        reporter = self._get_reporter()
        xmlpath = self._get_xml_path()
        compoundtree = ET.parse(xmlpath)
        root = compoundtree.getroot()
        if len(root) > 1:
            reporter.xml_assert(xmlpath, "more than one compound in a file")
        if root[0].tag != 'compounddef':
            reporter.xml_assert(xmlpath, "expected <compounddef> as the first tag")
            return
        briefelem = None
        detailselem = None
        missing_members = set(self._members.values())
        for elem in root[0]:
            if elem.tag == 'compoundname':
                if elem.text != self.get_name():
                    reporter.xml_assert(xmlpath,
                            "compound name mismatch: '{0}' (in index.xml) vs. '{1}'"
                            .format(self.get_name(), elem.text))
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
                    reporter.xml_assert(xmlpath,
                            "unknown inner compound type '{0}'".format(reftype))
            elif elem.tag == 'sectiondef':
                # TODO: Handle header and description elements
                kind = elem.attrib['kind']
                section = MemberSection(kind)
                self._sections.append(section)
                for memberelem in elem.iter('memberdef'):
                    refid = memberelem.attrib['id']
                    member = self._members[refid]
                    member.load_details_from_element(memberelem, xmlpath)
                    section.add_member(member)
                    if member in missing_members:
                        missing_members.remove(member)
                    # Enum values need special handling, but are not worth
                    # extra generalization.
                    if isinstance(member, Enum):
                        missing_members.difference_update(member.get_values())
            else:
                if not self._load_element(elem):
                    reporter.xml_assert(xmlpath,
                            "unknown compound child element '{0}'".format(elem.tag))
        if missing_members:
            reporter.xml_assert(xmlpath, 'members without section')
        self._process_descriptions(briefelem, detailselem, None)
        self._loaded = True

    def _unexpected_inner_compound(self, typename, compound):
        """Report a parsing error for an unexpected inner compound reference."""
        reporter = self._get_reporter()
        xmlpath = self._get_xml_path()
        reporter.xml_assert(xmlpath,
                "unexpected inner {0}: {1}".format(typename, compound))

    def _load_inner_file(self, compound):
        """Process a reference to an inner file.

        Derived classes should override the method if the compound type can
        contain files as nested compounds.
        """
        self._unexpected_inner_compound("file", compound)

    def _load_inner_dir(self, compound):
        """Process a reference to an inner directory.

        Derived classes should override the method if the compound type can
        contain directories as nested compounds.
        """
        self._unexpected_inner_compound("dir", compound)

    def _load_inner_group(self, compound):
        """Process a reference to an inner group.

        Derived classes should override the method if the compound type can
        contain groups as nested compounds.
        """
        self._unexpected_inner_compound("group", compound)

    def _load_inner_namespace(self, compound):
        """Process a reference to an inner namespace.

        Derived classes should override the method if the compound type can
        contain namespaces as nested compounds.
        """
        self._unexpected_inner_compound("namespace", compound)

    def _load_inner_class(self, compound):
        """Process a reference to an inner class.

        Derived classes should override the method if the compound type can
        contain classes as nested compounds.
        """
        self._unexpected_inner_compound("class", compound)

    def _load_element(self, element):
        """Load data from a child XML element.

        This method is called for all XML elements under the <compounddef>
        element that are not handled directly by the Compound class.
        Derived classes should return True if they process the element.
        """
        return False

    def get_groups(self):
        return self._groups

    def show_base(self):
        """Format information for common properties.

        This extends Entity.show_base() by adding properties that are common to
        all compounds.
        """
        Entity.show_base(self)
        if self._groups:
            print 'Groups:   {0}'.format(', '.join(map(str, self._groups)))

    def show_members(self):
        """Show list of members.

        This method is provided for use in show() methods of derived classes
        to print the list of members.
        """
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
        self._is_source_file = None

    def _load_inner_class(self, compound):
        compound.add_file(self)
        self._classes.add(compound)

    def _load_inner_namespace(self, compound):
        compound.add_file(self)
        self._namespaces.add(compound)

    def _load_element(self, elem):
        if elem.tag == 'location':
            self._path = elem.attrib['file']
            extension = os.path.splitext(self._path)[1]
            self._is_source_file = (extension in ('.c', '.cpp', '.cu'))
            return True
        return False

    def set_directory(self, directory):
        self._directory = directory

    def get_reporter_location(self):
        return reporter.Location(self._path, None)

    def get_path(self):
        return self._path

    def get_directory(self):
        return self._directory

    def is_source_file(self):
        return self._is_source_file

    def show(self):
        self.show_base()
        print 'Path:      {0}'.format(self._path)
        print 'Directory: {0}'.format(self._directory)
        print 'Source:    {0}'.format(self._is_source_file)
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

    def get_reporter_location(self):
        return reporter.Location(self._path, None)

    def get_path(self):
        return self._path

    def get_parent(self):
        return self._parent

    def get_subdirectories(self):
        return self._subdirs

    def show(self):
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

    def get_reporter_location(self):
        return self._doclocation.get_reporter_location()

    def is_anonymous(self):
        return 'anonymous_namespace{' in self.get_name()

    def show(self):
        self.show_base()
        print 'Doc. loc.: {0}'.format(self._doclocation.get_full_string())
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

    def get_reporter_location(self):
        return self._location.get_reporter_location()

    def get_files(self):
        return self._files

    def is_local(self):
        if len(self._files) > 1:
            return False
        for fileobj in self._files:
            if not fileobj.is_source_file():
                return False
        return True

    def show(self):
        self.show_base()
        print 'Namespace:  {0}'.format(self._namespace)
        if self._outerclass:
            print 'Outer cls:  {0}'.format(self._outerclass)
        location = self._location
        print 'Location:   {0}'.format(location.get_location().get_full_string())
        print 'Body loc:   {0}'.format(location.get_body_location().get_full_string())
        _show_list('Inner classes', self._innerclasses)
        self.show_members()

#####################################################################
# Top-level container class

def _get_compound_type_from_kind(kind):
    """Map compound kinds from Doxygen XML to internal class types."""
    if kind == 'file':
        return File
    elif kind == 'dir':
        return Directory
    elif kind == 'group':
        return Group
    elif kind == 'namespace':
        return Namespace
    elif kind in ('class', 'struct', 'union'):
        return Class
    else:
        return None

def _get_member_type_from_kind(kind):
    """Map member kinds from Doxygen XML to internal class types."""
    if kind == 'define':
        return Define
    elif kind == 'variable':
        return Variable
    elif kind == 'typedef':
        return Typedef
    elif kind == 'enum':
        return Enum
    elif kind == 'enumvalue':
        return EnumValue
    elif kind == 'function':
        return Function
    elif kind == 'friend':
        return FriendDeclaration
    else:
        return None

class DocumentationSet(object):

    """Root object for Doxygen XML documentation tree.

    On initialization, it reads the index.xml file from the Doxygen XML output,
    which contains the list of entities.  Only the ID and name for the entities,
    and the parent compounds for members, are available from this file.

    load_details() can be called to load the detailed compound XML files.
    This constructs relations between compound entities, and initializes other
    attributes for the entities.

    load_file_details() does the same as load_details(), except that it leaves
    those compound XML files unloaded that do not affect file objects or their
    parent hierarchy.  This saves some time if details for actual code
    constructs like namespaces, classes or members are not necessary.

    merge_duplicates() can then be called to remove members with different IDs,
    but that actually reference the same code entity.  For some reason, Doxygen
    seems to produce these in certain cases.
    """

    def __init__(self, xmlroot, reporter):
        """Initialize the documentation set and read index data."""
        self._xmlroot = xmlroot
        self._reporter = reporter
        xmlpath = os.path.join(xmlroot, 'index.xml')
        indextree = ET.parse(xmlpath)
        self._compounds = dict()
        self._members = dict()
        self._files = dict()
        for compoundelem in indextree.getroot():
            name = compoundelem.find('name').text
            refid = compoundelem.attrib['refid']
            kind = compoundelem.attrib['kind']
            if kind in ('page', 'example'):
                # TODO: Model these types as well
                continue
            compoundtype = _get_compound_type_from_kind(kind)
            if compoundtype is None:
                reporter.xml_assert(xmlpath,
                        "unknown compound kind '{0}'".format(kind))
                continue
            compound = compoundtype(name, refid)
            compound.set_documentation_set(self)
            self._compounds[refid] = compound
            for memberelem in compoundelem.iter('member'):
                name = memberelem.find('name').text
                refid = memberelem.attrib['refid']
                kind = memberelem.attrib['kind']
                if refid in self._members:
                    member = self._members[refid]
                    membertype = _get_member_type_from_kind(kind)
                    if not isinstance(member, membertype):
                        reporter.xml_assert(xmlpath,
                                "id '{0}' used for multiple kinds of members"
                                .format(refid))
                        continue
                else:
                    membertype = _get_member_type_from_kind(kind)
                    if membertype is None:
                        reporter.xml_assert(xmlpath,
                                "unknown member kind '{0}'".format(kind))
                        continue
                    member = membertype(name, refid)
                    member.set_documentation_set(self)
                    self._members[refid] = member
                member.add_parent_compound(compound)
                compound.add_member(member)

    def load_file_details(self, filelist=None):
        """Load detailed XML files for all files and possible parents of files.

        If filelist is set, it should be a list of file paths, and details will
        be loaded only for files in those paths.  The path format should match
        what Doxygen writes into the files (with Gromacs setup, it seems to be
        absolute paths)."""
        for compound in self._compounds.itervalues():
            if isinstance(compound, (Directory, Group)):
                compound.load_details()
            elif not filelist and isinstance(compound, File):
                compound.load_details()
                self._files[compound.get_path()] = compound
        if filelist:
            # We can't access the full path from the File object before the
            # details are loaded, because Doxygen does not write that into
            # index.xml.  But we can use the Directory objects (which were
            # loaded above) to get the path.
            for compound in self._compounds.itervalues():
                if isinstance(compound, File):
                    dirobj = compound.get_directory()
                    if not dirobj:
                        continue
                    abspath = compound.get_directory().get_path()
                    abspath = os.path.join(abspath, compound.get_name())
                    if abspath in filelist:
                        compound.load_details()
                        self._files[compound.get_path()] = compound

    def load_details(self):
        """Load detailed XML files for each compound."""
        for compound in self._compounds.itervalues():
            compound.load_details()
            if isinstance(compound, File):
                self._files[compound.get_path()] = compound
        # TODO: Add links to files using location

    def merge_duplicates(self):
        """Merge duplicate member definitions based on body location.

        At least for some functions that are declared in a header, but have
        their body in a source file, Doxygen seems to create two different IDs,
        but the contents of the members are the same, except for the location
        attribute.  This method merges members that have identical name and
        body location into a single member that keeps the information from both
        instances (they should only differ in the location attribute and in
        parent compounds).  Both IDs point to the merged member after this
        method.
        """
        members_by_body = dict()
        for member in self._members.itervalues():
            bodyloc = member.get_body_location()
            if bodyloc:
                index = (bodyloc, type(member), member.get_name())
                if index not in members_by_body:
                    members_by_body[index] = []
                members_by_body[index].append(member)
        for memberlist in members_by_body.itervalues():
            if len(memberlist) > 1:
                declaration = None
                otherdeclarations = []
                definition = None
                for member in memberlist:
                    if member.has_same_body_location():
                        if definition is not None:
                            self._reporter.xml_assert(None,
                                    "duplicate definition for a member '{0}'"
                                    .format(definition))
                            continue
                        definition = member
                    elif declaration is None:
                        declaration = member
                    else:
                        otherdeclarations.append(member)
                if otherdeclarations:
                    # TODO: gmx_cpuid.c produces some false positives
                    details = []
                    for otherdeclaration in otherdeclarations:
                        details.append('{0}: another declaration is here'
                                .format(otherdeclaration.get_reporter_location()))
                    details.append('{0}: definition is here'
                            .format(declaration.get_body_location()))
                    text = "duplicate declarations for a member '{0}'".format(declaration)
                    self._reporter.code_issue(declaration, text, details)
                    continue
                self._members[definition.get_id()] = declaration
                declaration.merge_definition(definition)
                for compound in definition.get_parent_compounds():
                    compound.replace_member(definition, declaration)

    def get_reporter(self):
        """Return reporter object to use for reporting issues.

        This method is used in the entity classes to access the reporter when
        they are parsing the XML files.
        """
        return self._reporter

    def get_xmlroot(self):
        """Return root of the Doxygen XML directory."""
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
        # self._members can contain duplicates because of merge_duplicates()
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

    def get_namespaces(self, name=None):
        if name:
            return self.get_compounds(Namespace, lambda x: x.get_name() in name)
        else:
            return self.get_compounds(Namespace)

    def get_classes(self, name=None):
        if name:
            return self.get_compounds(Class, lambda x: x.get_name() in name)
        else:
            return self.get_compounds(Class)

    def get_functions(self, name):
        return self.get_members(Member, lambda x: x.get_name() in name)

#####################################################################
# Code for running in script mode

def main():
    """Run the script in for debugging/Doxygen XML output inspection."""
    import sys

    from optparse import OptionParser

    from reporter import Reporter

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
    # TODO: Add option for other types, and make them work
    parser.add_option('-f', '--show-function', action='append',
                      help='Show details of given function')
    options, args = parser.parse_args()

    reporter = Reporter()

    sys.stderr.write('Loading index.xml...\n')
    docset = DocumentationSet(options.root_dir, reporter)
    reporter.write_pending()
    sys.stderr.write('Loading details...\n')
    docset.load_details()
    reporter.write_pending()
    sys.stderr.write('Processing...\n')
    docset.merge_duplicates()
    reporter.write_pending()

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
