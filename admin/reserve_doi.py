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

# This python script registers a doi for a source archive or manual
# of a GROMACS release version, to allow for proper referencing in
# the future.
# It is called from CMAKE and returns the final doi and the submission id
# CMAKE passes the type of registration to it,
# as well as the current version as the arguments


parser = argparse.ArgumentParser(description='Select doi type to register')
parser.add_argument('type',metavar='type', type=str,
        help='The type of doi that will be registered',
        choices=['manual', 'source'])
parser.add_argument('version',metavar='ver', type=str,
        help='The current GROMACS version')
parser.add_argument('secret',metavar='secret', type=str,
        help='File location for our secret text file on Jenkins, passed by CMake')


# get the kind of doi we will request
doi_type = (vars(parser.parse_args()))['type']
# get GROMACS version from input
version = (vars(parser.parse_args()))['version']
# get secret file path
secret_path = (vars(parser.parse_args()))['secret']

# set some general variables that are true for both cases
creators = [{'name' : 'Abraham, Mark', 'affiliation': 'KTH'}, {'name':'van der Spoel, David', 'affiliation':'Uppsala University, ICM'},{'name':'Hess, Berk','affiliation':'KTH'}, {'name':'Lindahl, Erik','affiliation':'KTH'}]
access_right = 'open'
license = 'LGPL-2.1'
# TODO add complete list of contributors
# contributors = [{}]
# TODO find a way to securely include the access token here
# read from file on Jenkins that is kept there in private
# currently for a local file on my computer for the sandbox :)
try:
    with open (secret_path, "r") as tokenfile:
        token="".join(line.rstrip() for line in tokenfile)

except IOError:
    print "Could not open file containing authentication token ",tokenfile
    exit(1)

zenodo_headers = {"Content-Type": "application/json"}
# now we can set all variables that depend on our doi type
# just one giant if statement

# TODO change this to a dictionary lookup, so it looks less ugly
if doi_type == "manual":
    upload_type = 'publication'
    publication_type = 'softwaredocumentation'
    title = 'GROMACS Reference Manual'
    description = 'GROMACS reference manual, version '+version
elif doi_type == "source":
    upload_type = 'software'
    publication_type = ''
    title = 'GROMACS source code'
    description = 'GROMACS source code tarball, version '+version

# create the actual object containing the information for the new release
# metadata for release
zenodo_metadata = {
        'metadata': {
            'title':title,
            'upload_type':upload_type,
            'creators':creators,
            'publication_type':publication_type,
            'access_right':access_right,
            'license':license,
            'title':title,
            'description':description,
            'version':version,
            'prereserve_doi':'true'
            }
        }
# generate the new object on zenodo
r = requests.post('https://zenodo.org/api/deposit/depositions',
        params={'access_token': token}, json={},
        headers=zenodo_headers)
if r.status_code != 201:
    print "Zenodo server returned incorrect status code"
    exit(1)

# get the id generated for this release
post_id = r.json()['id']
# set the metadata for the release
r = requests.put('https://zenodo.org/api/deposit/depositions/%s' % post_id,
        params={'access_token': token},headers=zenodo_headers,
        data=json.dumps(zenodo_metadata))
if r.status_code != 200:
    print "Zenodo server returned incorrect status code"
    exit(1)


# nested request to get to the doi in metadata->prereserve_doi
doi = ((r.json()['metadata'])['prereserve_doi'])['doi']

# return doi and id to be saved in CMake variables
print doi,post_id
