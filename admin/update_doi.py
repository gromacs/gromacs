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

import os
import argparse
import re
import requests
import json

# This python script request to update a resource on Zenodo
# CMake passes the type of registration to it,
# the software/manual version,
# as well as the id needed to update the resource
# the script returns the new doi and submission if for the new resource
# and deletes the old file found there


parser = argparse.ArgumentParser(description='Select doi type to register')
parser.add_argument('type',metavar='type', type=str,
        help='The type of doi that will be registered',
        choices=['manual', 'source'])
parser.add_argument('version',metavar='ver', type=str,
        help='The current GROMACS version')
parser.add_argument('id',metavar='id', type=int,
        help='The old submission id to be updated')
parser.add_argument('secret',metavar='secret', type=str,
        help='File location for our secret text file on Jenkins, passed by CMake')

# get the kind of doi we will request
doi_type = (vars(parser.parse_args()))['type']
# get GROMACS version from input
version = (vars(parser.parse_args()))['version']
# get the submission id from input
zenodo_id = (vars(parser.parse_args()))['id']
# get secret file path
secret_path = (vars(parser.parse_args()))['secret']

try:
    with open (secret_path, "r") as tokenfile:
        token="".join(line.rstrip() for line in tokenfile)

except IOError:
    print "Could not open file containing authentication token ",tokenfile
    exit(1)

# now we can set all variables that depend on our doi type

# TODO change this to a dictionary lookup, so it looks less ugly
if doi_type == "manual":
    upload_type = 'publication'
    publication_type = 'softwaredocumentation'
    title = 'Random Reference Manual'
    description = 'Random reference manual, version '+version
elif doi_type == "source":
    upload_type = 'software'
    publication_type = ''
    title = 'Random source code'
    description = 'Random source code tarball, version '+version

zenodo_headers = {"Content-Type": "application/json"}
# file to be uploaded
zenodo_data = {'filename':file_path}

# update the old information to get new version
r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions/%s/actions/newversion' % zenodo_id,
        params={'access_token': token})

if r.status_code != 201:
    print "Zenodo server returned incorrect status code"
    exit(1)

# get url for the new draft, needed to get new doi and id
draft_url = (r.json()['links'])['latest_draft']

# generate new requests object with updated information
# but keep the old one
r_update = requests.get(draft_url,params={'access_token': token})

if r.status_code != 200:
    print "Zenodo server returned incorrect status code"
    exit(1)

# get new zenodo id to upload the files to
zenodo_id = r_update.json()['id']

# get old metadata to update
zenodo_metadata = r_update.json()['metadata']

# TODO should also keep option to update other fields here
# update entry metadata
zenodo_metadata['version'] = version
zenodo_metadata['description'] = description

# update full entry data with new metadata
new_data = r_update.json()
new_data['metadata'] = zenodo_metadata

# upload new data
r_update = requests.put('https://sandbox.zenodo.org/api/deposit/depositions/%s' % zenodo_id,
        params={'access_token': token},headers=zenodo_headers,
        data=json.dumps(new_data))

if r.status_code != 200:
    print "Zenodo server returned incorrect status code"
    exit(1)

# nested request to get to the doi in metadata->prereserve_doi
doi = ((r_update.json()['metadata'])['prereserve_doi'])['doi']

# remove old files to make it possible to upload the new ones
# get file object
r_update = requests.get('https://sandbox.zenodo.org/api/deposit/depositions/%s/files' % zenodo_id,
     params={'access_token': token})

if r.status_code != 200:
    print "Zenodo server returned incorrect status code"
    exit(1)

# file_list is a list of attributes, need to get id value from it
file_list = r_update.json()
# get first element and extract file id from it
file_id = (file_list[0])['id']
# actually delete the old file
r_update = requests.delete('https://sandbox.zenodo.org/api/deposit/depositions/%s/files/%s' % (zenodo_id, file_id),
     params={'access_token': token})

if r.status_code != 204:
    print "Zenodo server returned incorrect status code"
    exit(1)

# return doi and id to be saved in CMake variables
print doi,zenodo_id




