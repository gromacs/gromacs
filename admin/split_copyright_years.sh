#!/usr/bin/env bash
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# Finds copyright statements that have more than five years, such as
#
#   ... Copyright (c) 2012,2013,2016,2017,2018,2019, by the GROMACS ...
#
# and splits them into multiple lines, like
#
#   ... Copyright (c) 2012,2013,2016,2017,2018, by the GROMACS development team.
#   ... Copyright (c) 2019, by the GROMACS development team, led by
#
# so that the copyright checker recognizes the second line as one that can be extended.
# This will need to be re-rerun every few years as the length of the last line grows.

sed -i 's/\( \?.\) Copyright (c) \([0-9][0-9][0-9][0-9],[0-9][0-9][0-9][0-9],[0-9][0-9][0-9][0-9],[0-9][0-9][0-9][0-9],[0-9][0-9][0-9][0-9]\),\([0-9][0-9][0-9][0-9],.*,\) by the GROMACS/\1 Copyright (c) \2 by the GROMACS development team.\n\1 Copyright (c) \3 by the GROMACS/g' $(git grep -l "by the GROMACS")
