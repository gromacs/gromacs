#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
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

"""Script for generating spherical test configurations."""

import math
import random
import sys

def dot(a, b):
    return sum(x*y for (x, y) in zip(a, b))

def norm(a):
    return math.sqrt(dot(a, a))

def angle(a, b):
    return math.degrees(math.acos(dot(a, b) / (norm(a) * norm(b))))

def minangle(a, list):
    minangle = 180
    minindex = -1
    for index, x in enumerate(list):
        xangle = angle(a, x)
        if xangle < minangle:
            minangle = xangle
            minindex = index
    return (minangle, minindex)

def get_single_vec():
    while True:
        x = random.randint(-1000, 1000) / 1000.0
        y = random.randint(-1000, 1000) / 1000.0
        z = random.randint(-1000, 1000) / 1000.0
        pos = (x, y, z)
        dist = norm(pos)
        if dist <= 1.0 and dist > 0.25:
            return pos

def write_gro_title(fp, title, atomcount):
    fp.write(title + '\n')
    fp.write('%5d\n' % (atomcount))

def write_gro_atom(fp, resnr, resname, atomname, index, x):
    fp.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %
            (resnr, resname, atomname, index, x[0], x[1], x[2]))

def write_gro_box(fp, box):
    fp.write('%10.5f%10.5f%10.5f\n' % box)

random.seed(1097)
center = (0, 0, 0)
cutoff = 20
possamples = 500
negsamples = 500
totsamples = 10000

sys.stderr.write("Generating reference points\n")
refpoints = []
refpoints.append((0, 0, -1))
refpoints.append((-0.5, 0.6, 0.1))
refpoints.append((-0.5, -0.5, 0.25))
while len(refpoints) < 30:
    pos = get_single_vec()
    if pos[0] > 0 and pos[1] > pos[0] and pos[2] > 0:
        refpoints.append(pos)

sys.stderr.write("Generating test points\n")
postestpoints = []
negtestpoints = []
hits = 0
samplecount = 0
while samplecount < totsamples or len(postestpoints) < possamples or len(negtestpoints) < negsamples:
    pos = get_single_vec()
    (pangle, index) = minangle(pos, refpoints)
    if pangle < cutoff:
        hits += 1
        if len(postestpoints) < possamples:
            postestpoints.append(pos)
    if pangle > cutoff:
        if len(negtestpoints) < negsamples:
            negtestpoints.append(pos)
    samplecount += 1

cfrac = float(hits) / samplecount
errest = math.sqrt((cfrac - cfrac*cfrac) / samplecount)
sys.stderr.write('Cutoff: %f angles\n' % (cutoff))
sys.stderr.write('Estimated covered fraction: %f +- %f\n' % (cfrac, errest))

debugfp = open('debug.txt', 'w')
fp = sys.stdout
count = 1 + len(refpoints) + len(postestpoints) + len(negtestpoints)
write_gro_title(fp, 'Spherical test case, cutoff %f, cfrac %f +- %f' %
        (cutoff, cfrac, errest) , count)
n = 1
write_gro_atom(fp, 1, 'C', 'C', n, center)
n += 1
for i in range(len(refpoints)):
    write_gro_atom(fp, 2, 'R', 'R', n, refpoints[i])
    n += 1
for i in range(len(postestpoints)):
    write_gro_atom(fp, 3, 'TP', 'TP', n, postestpoints[i])
    x = postestpoints[i]
    xangle, index = minangle(x, refpoints)
    refx = refpoints[index]
    debugfp.write('%3d%8.3f%8.3f%8.3f  %4.1f  %2d%8.3f%8.3f%8.3f\n' %
            (n-1, x[0], x[1], x[2], xangle, index, refx[0], refx[1], refx[2]))
    n += 1
for i in range(len(negtestpoints)):
    write_gro_atom(fp, 4, 'TN', 'TN', n, negtestpoints[i])
    n += 1
write_gro_box(fp, (10, 10, 10))
debugfp.close()
