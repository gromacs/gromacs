#!/usr/bin/env python3
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2022- The GROMACS Authors
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

"""Script to dispatch a post-merge-acceptance pipeline.

Submits a new post-merge pipeline for the specified branch, useful for
testing issues that only arise in those jobs.

Requires the `python_gitlab <https://python-gitlab.readthedocs.io/en/stable/index.html>`__
package for Python access to the GitLab API.

Author:
    * Paul Bauer <paul.bauer.q@gmail.com>
"""

import argparse

import gitlab


def submit_gitlab_pipeline(auth_token, branch='main'):
    """Submit a post merge pipeline to GitLab server.

    The pipeline is by default created for main, but can be run for any branch.

    Needs access to gitlab python API. For this, someone running the script needs
    to have installed the python library (pip install python-gitlab) and created
    a personal access token that grants access to the GROMACS project
    (https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html).

    Throws if the authentication process, or the submission of the pipeline fails.
    """
    gl = gitlab.Gitlab('https://gitlab.com', private_token=auth_token)
    # The project ID for GROMACS is hardcoded here
    project_id = 17679574
    project = gl.projects.get(project_id)

    print('Going to start pipeline with following arguments')
    print(r'Branch = ', branch)
    print(r'POST_MERGE_ACCEPTANCE = true')

    return project.pipelines.create({
        'ref': branch,
        'variables': [{'key': 'POST_MERGE_ACCEPTANCE', 'value': 'true'}]})


parser = argparse.ArgumentParser(
    description='Options for submitting post-merge pipelines.')

parser.add_argument('--token', type=str, required=True,
                    help='GitLab access token needed to launch pipelines')

parser.add_argument('--branch', type=str, default='main',
                    help='Branch to run pipeline for (default "main")')

if __name__ == "__main__":
    args = parser.parse_args()
    submit_gitlab_pipeline(args.token, args.branch)
