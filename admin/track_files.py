#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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

import datetime
import subprocess
import sys

from optparse import OptionParser
# from subprocess import Popen, PIPE

def convert_date(git_date_string):
    year, month, day = git_date_string.split('-')
    return datetime.date(int(year), int(month), int(day))

class CommitInfo(object):
    def __init__(self, sha1, commit_date, author_date, subject):
        self._hash = sha1
        self._commit_date = convert_date(commit_date)
        self._author_date = convert_date(author_date)
        self._subject = subject

    def get_commit_year(self):
        return self._commit_date.year

    def get_oneline(self):
        return '{0} {1:.50}'.format(self._hash, self._subject)

class FileInfo(object):
    def __init__(self):
        self._commits = []
        self._creation_commit = None
        self._years = set()

    def set_creation_commit(self, commit):
        self._creation_commit = commit
        self.add_commit(commit)

    def add_commit(self, commit):
        self._commits.append(commit)
        self._years.add(commit.get_commit_year())

    def get_years(self):
        return ','.join(sorted([str(x) for x in self._years]))

    def get_creation_commit(self):
        return self._creation_commit

    def get_commits(self):
        return self._commits

def parse_changes(lines, current_files, added_files, modified_files):
    while lines and not lines[0]:
        lines.pop(0)
    while lines and lines[0].startswith(':'):
        parts = lines.pop(0).split('\t')
        status = parts[0].split(' ')[4]
        filename = parts[1]
        if status == 'A':
            if filename in current_files:
                added_files.add(current_files.pop(filename))
        elif status == 'M':
            if filename in current_files:
                modified_files.add(current_files[filename])
        elif status.startswith('R'):
            srcname = filename
            destname = parts[2]
            if destname in current_files:
                fileinfo = current_files.pop(destname)
                modified_files.add(fileinfo)
                current_files[srcname] = fileinfo

def parse_commit(lines, current_branches, files, excludelist):
    (sha1, parents, commit_date, author_date, subject) = lines.pop(0).split(';', 4)
    parents = parents.split()
    commit = CommitInfo(sha1, commit_date, author_date, subject)
    if not current_branches:
        current_branches[sha1] = files.copy()
    elif not sha1 in current_branches:
        sys.stderr.write('Unknown commit: {0} {1}\n'.format(sha1, subject))
        sys.exit(1)
    current_files = current_branches.pop(sha1)
    if len(parents) > 1:
        branch_files = current_files.copy()
    added_files = set()
    modified_files = set()
    parse_changes(lines, current_files, added_files, modified_files)
    current_branches[parents[0]] = current_files
    parents.pop(0)
    for parent in parents:
        if not lines or not lines[0].startswith(sha1):
            break
        lines.pop(0)
        current_files = branch_files.copy()
        current_branches[parent] = current_files
        parse_changes(lines, current_files, added_files, modified_files)
    for exclude in excludelist:
        if exclude.startswith(sha1):
            return
    added_files.difference_update(modified_files)
    for added in added_files:
        added.set_creation_commit(commit)
    for modified in modified_files:
        modified.add_commit(commit)

def main():
    parser = OptionParser()
    parser.add_option('--stdin', action='store_true',
                      help='Read list of files from stdin')
    parser.add_option('--rev', default='HEAD',
                      help='Set the revision to track back from')
    parser.add_option('--since',
                      help='Set the date to track back to')
    parser.add_option('--exclude-commits',
                      help='Provide a file with a list of commits to ignore')
    parser.add_option('--track-merges', action='store_true',
                      help='Track through merges')
    parser.add_option('--created', action='store_true',
                      help='Report when the files were created')
    parser.add_option('--modified', action='store_true',
                      help='Report when the files were modified')
    parser.add_option('--commits', action='store_true',
                      help='Report commit information')
    options, args = parser.parse_args()

    if options.stdin:
        filelist = sys.stdin.read().splitlines()
    else:
        lscmd = ['git', 'ls-tree', '--name-only', '-r', options.rev, '--']
        lscmd.extend(args)
        filelist = subprocess.check_output(lscmd).splitlines()
    excludelist = []
    if options.exclude_commits:
        with open(options.exclude_commits, 'r') as excludefile:
            excludelist = excludefile.read().splitlines()
    files = dict()
    for filename in filelist:
        files[filename] = FileInfo()
    logcmd = ['git', 'log', '-m', '-M', '--raw',
              '--since=' + options.since, '--date=short', '--format=%h;%p;%cd;%ad;%s']
    if not options.track_merges:
        logcmd.append('--first-parent')
    logcmd.append(options.rev)
    logoutput = subprocess.check_output(logcmd)
    lines = logoutput.splitlines()
    current_branches = dict()
    while lines:
        parse_commit(lines, current_branches, files, excludelist)

    for filename in filelist:
        fileinfo = files[filename]
        summary = []
        details = []
        if options.created:
            creation_commit = fileinfo.get_creation_commit()
            if creation_commit:
                summary.append('created {0}'.format(creation_commit.get_commit_year()))
                if options.commits:
                    details.append('created: {0}'.format(creation_commit.get_oneline()))
        if options.modified:
            years = fileinfo.get_years()
            commits = fileinfo.get_commits()
            if commits:
                summary.append('modified {0}'.format(years))
                if options.commits:
                    for commit in commits:
                        details.append('{0}: {1}'.format(commit.get_commit_year(), commit.get_oneline()))
        if summary:
            sys.stdout.write('{0}: {1}\n'.format(filename, ', '.join(summary)))
            if details:
                sys.stdout.write('  ' + '\n  '.join(details) + '\n')

if __name__ == "__main__":
    main()
