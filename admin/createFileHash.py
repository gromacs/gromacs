#! /usr/bin/env python
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
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
import hashlib, hmac, os, stat, sys, re
from re import search

"""
Calculate hash of files in build tree to allow checking against
stored hashes in case of the tree not being in git (e.g. if the
program is build from a release tarball.

Based on example script found here:
    https://unix.stackexchange.com/a/35847
"""

def is_in_whitelist(name):
    """Return true if file is white listed to be included in hash calculation."""
    in_whitelist = False
    whitelist = ["\.cpp$", "\.h$", "\.cuh$", "\.cu$", "\.clh$", "CMakeList.txt$", "\.cmake$", "\.in$", "\.cmakein$", "\.py$"]
    for item in whitelist:
        if search(item, name):
            in_whitelist = True
            break

    return in_whitelist

def is_blacklisted(name):
    """Return if a file has been explicitly blacklisted.

    """
    is_blacklisted = False
    blacklist = ["gmx-completion"]
    for item in blacklist:
        if search(item, name):
            is_blacklisted = True
            break

    return is_blacklisted

def file_hash(name):
    """Return the hash of the contents of the specified file, as a hex string

    Reads file in chunks of 16384 bytes and calculates the hash of the complete
    file afterwards.
    The hashing algorithm used is sha256, to avoid accidental clashes when using
    a more simple algorithm such as md5.
    """
    f = open(name, 'rb')
    h = hashlib.sha256()
    while True:
        buf = f.read(16384)
        if len(buf) == 0: break
        h.update(buf)
    f.close()
    return h.hexdigest()

def traverse(h, path, original_path):
    """Recursive function to traverse a file path until a regular file is found.
    Walks down the path given as the input and updates the hash function with
    information of new files that are found on bottom of the list.

    Information used to calculate the hash are the name and the contents of the file.
    Uses both absolute and relative path to make sure only the relative path is used
    to calculate the hash.

    Ignores files that are not in the white-list and also skips files that are
    explicitly blacklisted.
    Other things that are ignored are symlinks and all kinds of special files.
    """
    rs = os.lstat(path)
    quoted_name = repr(os.path.relpath(path, original_path))
    if stat.S_ISDIR(rs.st_mode):
        for entry in sorted(os.listdir(path)):
            traverse(h, os.path.join(path, entry), original_path)
    elif stat.S_ISREG(rs.st_mode):
        # Only test files that actually take part in building GROMACS
        if (is_in_whitelist(path) and not is_blacklisted(path)):
            fullname = 'reg ' + quoted_name + ' '
            fullname += str(rs.st_size) + ' '
            fullname += file_hash(path) + '\n'
            h.update(fullname.encode('utf-8'))
    else: pass # silently symlinks and other special files

def main():
    """Run the hashing script.

    Takes single directory to hash files in.

    """
    import os
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Hash all white listed files in a single directory')
    parser.add_argument('-s',
            '--source-root',
            help='Source tree directory, can be specified multiple times to get several directories hashed',
            nargs='*',
            required=True)
    parser.add_argument('-o',
            '--output-file',
            help='File to write hash to.',
            default='hashresult')

    args = parser.parse_args()
    
    outfile_path = args.output_file
    h = hashlib.sha256()
    for input_sources in args.source_root:
        traverse(h, input_sources, input_sources)
    
    end = 'end\n'
    h.update(end.encode('utf-8'))
    outputfile = open(outfile_path, 'w')
    outputfile.write(h.hexdigest())

if __name__ == '__main__':
    main()

