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

"""Script to dispatch and run release pipeline.

When run locally, submits a build pipeline to the
GitLab server with the necessary information to later
automatically upload generated artefacts to the relevant servers.

The script has two modes, local mode submits the job pipeline, while
in server mode we sync artefacts with the server.
Requires the `python_gitlab <https://python-gitlab.readthedocs.io/en/stable/index.html>`__
package for Python access to the GitLab API.

Author:
    * Andrey Alekseenko <al42and@gmail.com>
    * Mark Abraham <mark.j.abraham@gmail.com>
    * Paul Bauer <paul.bauer.q@gmail.com>
    * Eric Irrgang <ericirrgang@gmail.com>
"""

import argparse
from pathlib import Path
import os
import re
from shutil import copyfile
import subprocess
import typing


def submit_gitlab_pipeline(auth_token, ssh_key, branch='main', release_build=False, dry_run=True):
    """Submit a pipeline to GitLab server to run a release build.

    Authenticates user using token to the project, and tries to create new GROMACS_RELEASE
    pipeline. If release_build is true, and dry_run is false, this will generate a pipeline
    that will upload real artefacts to the manual and ftp servers if run from a release branch
    (or main).  Otherwise, the generated pipeline will always only upload files to a
    test location.

    The pipeline is by default created for main, but can be run for any branch.

    Needs access to gitlab python API. For this, someone running the script needs
    to have installed the python library (pip install python-gitlab) and created
    a personal access token that grants access to the GROMACS project
    (https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html).

    Throws if the authentication process, or the submission of the pipeline fails.
    """
    import gitlab

    gl = gitlab.Gitlab('https://gitlab.com', private_token=auth_token)
    # The project ID for GROMACS is hardcoded here
    project_id = 17679574
    project = gl.projects.get(project_id)

    # Set some variables
    release_build_str = f"{release_build}".lower()
    dry_run_str = f"{dry_run}".lower()

    ssh_key_string = Path(ssh_key).read_text()

    # add trailing new line for the ssh-key
    if ssh_key_string[:-1] != '\n':
        ssh_key_string += '\n'

    print('Going to start pipeline with following arguments')
    print(r'Branch = ', branch)
    print(r'GROMACS_RELEASE = ', release_build_str)
    print(r'DRY_RUN = ', dry_run_str)
    #print(r'SSH_PRIVATE_KEY = ', ssh_key_string)

    pipeline = project.pipelines.create({'ref': branch, 'variables': [{'key': 'GROMACS_RELEASE', 'value': release_build_str}, {
                                        'key': 'SSH_PRIVATE_KEY', 'value': ssh_key_string}, {'key': 'DRY_RUN', 'value': dry_run_str}]})


def upload_files(*, path: typing.Union[str, Path], options: typing.List[str], server: str, dry_run=False):
    """Actual upload command.

    Takes care of moving files to the server for public consumption.
    Needs to have successfully set up SSH key infrastructure before.
    """

    upload_command = ['rsync', '-rlptvP']
    upload_command.extend(options)

    upload_command.append(str(path))

    upload_command.append(server)

    print(upload_command)

    if not dry_run:
        ret = subprocess.run(upload_command, capture_output=True)
        if ret.returncode != 0:
            print(ret.stdout)
            print(ret.stderr)
            exit(ret.returncode)


def upload_release_artifacts():
    """Upload files to server.

    Sets up infrastructure for ssh to be able to authenticate to server and upload files.
    Key and host information is read from the GitLab environment variables that hold them
    after submitting the job.

    The following need to be true for uploading files to the actual location on
    the GROMACS project ftp and webservers:
        * Pipeline is run for either main, or one of the release branches
        * The pipeline is invoked with "GROMACS_RELEASE" set to "true"
        * The "DRY_RUN" variable needs to be explicitly set to "false"

    If any of those conditions above is not met, files are instead uploaded to a test
    location on the server (path prefixed with ".ci-test"), or, if "DRY_RUN" is not "false",
    no upload is attempted.

    For uploading a new version, the script will never overwrite existing files, as
    an additional measure of ensuing no problems with existing, world visible, files.
    """

    # Please see the variable documentation in docs/dev-manual/gitlab-ci.rst for explanation
    # of GROMACS specific  variables.
    dry_run_str = os.getenv('DRY_RUN')
    release_build_str = os.getenv('GROMACS_RELEASE')
    branch = os.getenv('CI_COMMIT_BRANCH')
    build_dir = os.getenv('BUILD_DIR')
    version = os.getenv('VERSION')

    if dry_run_str.lower() not in ('true', 'false'):
        raise ValueError(
            f'Wrong value of DRY_RUN: "{dry_run_str}". Only "true" and "false" are allowed')

    # we keep track of whether any command has failed so far
    is_upload = False
    dry_run = True
    if dry_run_str == 'false':
        dry_run = False

    full_version = version
    if not release_build_str == 'true':
        full_version += '-dev'

    overwrite_str = ''
    upload_location = './.ci-test'
    if re.match(r'^release-\d{4}$', branch) or re.match(r'^main$', branch):
        if release_build_str == 'true':
            is_upload = True
            overwrite_str = '--ignore-existing'
            upload_location = '.'
            # Only upload to real location if all preconditions are set.

    current_dir = os.getcwd()

    # set up manual front page repo
    os.mkdir('manual-front-page')
    os.chdir('manual-front-page')
    ret_init = subprocess.run(['git', 'init'], capture_output=True)
    ret_fetch = subprocess.run(
        ['git', 'fetch', 'https://gitlab.com/gromacs/deployment/manual-front-page.git', 'main'], capture_output=True)
    ret_checkout = subprocess.run(
        ['git', 'checkout', '-qf', 'FETCH_HEAD'], capture_output=True)
    os.chdir(current_dir)
    if ret_init.returncode != 0:
        print(ret_init.stdout)
        print(ret_init.stderr)
        exit(1)
    if ret_fetch.returncode != 0:
        print(ret_fetch.stdout)
        print(ret_fetch.stderr)
        exit(1)
    if ret_checkout.returncode != 0:
        print(ret_checkout.stdout)
        print(ret_checkout.stderr)
        exit(1)

    ret = subprocess.run(['make', 'html'], capture_output=True,
                         cwd=Path.cwd()/'manual-front-page')

    if ret.returncode != 0:
        print(ret.stdout)
        print(ret.stderr)
        exit(1)

    ftp_server = 'ftpbot@ftp.gromacs.org'
    manual_server = 'manualbot@manual.gromacs.org'

    website_file_location = build_dir+'/docs/html'
    frontpage_file_location = 'manual-front-page/_build/html'
    manual_file = f'manual-{version}.pdf'
    manual_file_location = website_file_location+'/'+manual_file

    website_path_on_server = f'{manual_server}:{upload_location}/{full_version}/'
    source_path_on_server = f'{ftp_server}:{upload_location}/gromacs/'
    regressiontests_path_on_server = f'{ftp_server}:{upload_location}/regressiontests/'
    manual_path_on_server = f'{ftp_server}:{upload_location}/manual/'
    frontpage_path_on_server = f'{manual_server}:{upload_location}/'

    website_upload_options = ['--chmod=u+rwX,g+rwX,o+rX']
    file_upload_options = ['--chmod=u+rw,g+rw,o+r']

    frontpage_upload_options = website_upload_options + \
        ['--exclude', '_sources', '--exclude',
            '.buildinfo', '--exclude', 'objects.inv']

    if is_upload and overwrite_str:
        website_upload_options.append(overwrite_str)
        file_upload_options.append(overwrite_str)

    manual_file = f'manual-{full_version}.pdf'
    copyfile(manual_file_location, manual_file)
    source_tarball = f'gromacs-{full_version}.tar.gz'
    regressiontests_tarball = f'regressiontests-{full_version}.tar.gz'

    os.chdir(website_file_location)
    upload_files(path='./',
                 options=website_upload_options,
                 server=website_path_on_server,
                 dry_run=dry_run)
    os.chdir(current_dir)
    upload_files(path=manual_file,
                 options=file_upload_options,
                 server=manual_path_on_server,
                 dry_run=dry_run)
    upload_files(path=source_tarball,
                 options=file_upload_options,
                 server=source_path_on_server,
                 dry_run=dry_run)
    upload_files(path=regressiontests_tarball,
                 options=file_upload_options,
                 server=regressiontests_path_on_server,
                 dry_run=dry_run)
    os.chdir(frontpage_file_location)
    upload_files(path='./',
                 options=frontpage_upload_options,
                 server=frontpage_path_on_server,
                 dry_run=dry_run)
    os.chdir(current_dir)


parser = argparse.ArgumentParser(
    description='Automatic release options.')

mode_group = parser.add_mutually_exclusive_group(required=True)
mode_group.add_argument('--local', action='store_true',
                        help='Set when running in local (submit pipeline) mode.')
mode_group.add_argument('--server', action='store_const', const=False, dest='local',
                        help='Set when running in server (upload artefacts) mode.')

parser.add_argument('--token', type=str,
                    help='GitLab access token needed to launch pipelines')
parser.add_argument('--ssh-key', type=str,
                    help='Path to SSH key needed to upload things to server. Pass in local mode to have it during the job')

release_parser = parser.add_mutually_exclusive_group(required=False)
release_parser.add_argument('--release', dest='release', action='store_true')
release_parser.add_argument(
    '--no-release', dest='release', action='store_const', const=False)

dry_run_parser = parser.add_mutually_exclusive_group(required=False)
dry_run_parser.add_argument('--dry-run', dest='dry_run', action='store_true')
dry_run_parser.add_argument(
    '--no-dry-run', dest='dry_run', action='store_const', const=False)

parser.add_argument('--branch', type=str, default='main',
                    help='Branch to run pipeline for (default "main")')


if __name__ == "__main__":

    args = parser.parse_args()

    if args.local:
        if args.token is None or args.branch is None or args.release is None or args.dry_run is None or args.ssh_key is None:
            raise RuntimeError(
                'Need to provide all command line options for running in local mode')

        ssh_key_location = Path(args.ssh_key).resolve()

        submit_gitlab_pipeline(args.token, ssh_key_location,
                               args.branch, args.release, args.dry_run)

    else:
        # --server mode
        upload_release_artifacts()
