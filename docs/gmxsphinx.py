#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

from sphinx import addnodes

class MdpNodeParser(object):
    def __init__(self):
        self._current_option = None

    def parse_option(self, env, text, nodes):
        nodes += addnodes.desc_name(text, text)
        self._current_option = text
        return text

    def parse_value(self, env, text, nodes):
        nodes += addnodes.desc_name(text, text)
        if self._current_option is None:
            return text
        return self._current_option + '=' + text

def setup(app):
    mdp_parser = MdpNodeParser()
    app.add_object_type('mdp', 'mdp',
            indextemplate='pair: %s; mdp option',
            parse_node = mdp_parser.parse_option,
            objname='mdp option')
    app.add_object_type('mdp-value', 'mdp-value',
            parse_node = mdp_parser.parse_value,
            objname='mdp value')
    app.add_object_type('cmake', 'cmake',
            indextemplate='pair: %s; cmake option',
            objname='CMake cache variable')
