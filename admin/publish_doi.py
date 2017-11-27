#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

try:
    import os
except ImportError:
    print "Could not import os"
    exit(1)
try:
    import argparse
except ImportError:
    print "Could not import argparse"
    exit(1)
try:
    import re
except ImportError:
    print "Could not import re"
    exit(1)
try:
    import requests
except ImportError:
    print "Could not import requests"
    exit(1)
try:
    import json
except ImportError:
    print "Could not import json"
    exit(1)

# This python script publishes the source or manual on Zenodo
# CMake passes the type of registration to it,
# the software/manual version,
# as well as the id needed to assign the documents


parser = argparse.ArgumentParser(description='Select doi type, set id, GROMACS version and path to files')
parser.add_argument('type',metavar='type', type=str,
        help='The type of material to be uploaded',
        choices=['manual', 'source'])
parser.add_argument('version',metavar='ver', type=str,
        help='The current GROMACS version')
parser.add_argument('id',metavar='id', type=int,
        help='The submission id')
parser.add_argument('path',metavar='path', type=str,
        help='Full path to current build directory')

# get the kind of doi we will request
doi_type = (vars(parser.parse_args()))['type']
# get GROMACS version from input
version = (vars(parser.parse_args()))['version']
# get the submission id from input
zenodo_id = (vars(parser.parse_args()))['id']
# get file path to this build from input
file_path = (vars(parser.parse_args()))['path']
# get secret file path from environment variable
secret_path = os.path.expandvars('$ZenodoTokenFile')

try:
    with open (secret_path, "r") as tokenfile:
        token="".join(line.rstrip() for line in tokenfile)

except IOError:
    print "Could not open file containing authentication token ",tokenfile
    exit(1)

# now we can set all variables that depend on our doi type
# just one if statement

if doi_type == 'manual':
        filename = 'manual-'+version+'.pdf'
        file_path = file_path+'/docs/html/'+filename
elif doi_type == 'source':
        filename = 'gromacs-'+version+'.tar.gz'
        file_path = file_path+'/'+filename

zenodo_files = ""
# those are static again and depend on the choices above
# this is just here to test that the file exists
# TODO test return value from the script and issue CMake FATAL_ERROR
# if something goes wrong
try:
    zenodo_files = {'file':open(file_path, 'rb')}
except IOError:
    print "Could not open file ",file_path
    exit(1)

zenodo_headers = {"Content-Type": "application/json"}
# file to be uploaded
zenodo_data = {'filename':file_path}

r = requests.get('https://zenodo.org/api/deposit/depositions/%s' % zenodo_id,
        params={'access_token': token})

bucket_url = r.json()['links']['bucket']

# upload file as pointed out here
# https://github.com/zenodo/zenodo/issues/833#issuecomment-324760423
r = requests.put('%s/%s' % (bucket_url,file_path),
        data=open(file_path, 'rb'),
        headers={"Accept":"application/json",
        "Authorization":"Bearer %s" % token,
        "Content-Type":"application/octet-stream"})

if r.status_code != 200:
    print "Zenodo server returned incorrect status code"
    exit(1)

# currently dummied out before we go live

# publish the material!
# r = requests.post('https://zenodo.org/api/deposit/depositions/%s/actions/publish' % zenodo_id,
#         params={'access_token': token})

# if r.status_code != 202:
#     print "Zenodo server returned incorrect status code"
#     exit(1)


