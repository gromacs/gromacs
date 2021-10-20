#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020,2021, by the GROMACS development team, led by
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

"""A utility module to help manage the matrix of configurations for CI testing and build containers.

When called as a stand alone script, prints a Docker image name based on the
command line arguments. The Docker image name is of the form used in the GROMACS
CI pipeline jobs.

Example::

    $ python3 -m utility --llvm --doxygen
    gromacs/ci-ubuntu-20.04-llvm-7-docs

See Also:
    :file:`buildall.sh`

As a module, provides importable argument parser and docker image name generator.

Note that the parser is created with ``add_help=False`` to make it friendly as a
parent parser, but this means that you must derive a new parser from it if you
want to see the full generated command line help.

Example::

    import utility.parser
    # utility.parser does not support `-h` or `--help`
    parser = argparse.ArgumentParser(
        description='GROMACS CI image creation script',
        parents=[utility.parser])
    # ArgumentParser(add_help=True) is default, so parser supports `-h` and `--help`

See Also:
    :file:`scripted_gmx_docker_builds.py`

Authors:
    * Paul Bauer <paul.bauer.q@gmail.com>
    * Eric Irrgang <ericirrgang@gmail.com>
    * Joe Jordan <e.jjordan12@gmail.com>
    * Mark Abraham <mark.j.abraham@gmail.com>
    * Gaurav Garg <gaugarg@nvidia.com>

"""

import argparse


parser = argparse.ArgumentParser(description='GROMACS CI image slug options.', add_help=False)
"""A parent parser for tools referencing image parameters.

This argparse parser is defined for convenience and may be used to partially initialize
parsers for tools.

.. warning:: Do not modify this parser.

    Instead, inherit from it with the *parents* argument to :py:class:`argparse.ArgumentParser`
"""

parser.add_argument('--cmake', nargs='*', type=str, default=['3.16.3', '3.17.2', '3.18.4', '3.21.2'], # new minimum required versions
                    help='Selection of CMake version to provide to base image')

compiler_group = parser.add_mutually_exclusive_group()
compiler_group.add_argument('--gcc', type=int, nargs='?', const=7, default=7,
                            help='Select GNU compiler tool chain. (Default) '
                                 'Some checking is implemented to avoid incompatible combinations')
compiler_group.add_argument('--llvm', type=str, nargs='?', const='7', default=None,
                            help='Select LLVM compiler tool chain. '
                                 'Some checking is implemented to avoid incompatible combinations')
# TODO currently the installation merely gets the latest beta version of oneAPI,
# not a specific version. GROMACS probably doesn't need to address that until
# oneAPI makes an official release. Also, the resulting container is a mix
# of packages with different betaXY version numbers, which hopefully works and
# is what Intel intends...
compiler_group.add_argument('--oneapi', type=str, nargs='?', const="2021.1.1", default=None,
                            help='Select Intel oneAPI package version.')

linux_group = parser.add_mutually_exclusive_group()
linux_group.add_argument('--ubuntu', type=str, nargs='?', const='20.04', default='20.04',
                         help='Select Ubuntu Linux base image. (default: ubuntu 20.04)')
linux_group.add_argument('--centos', type=str, nargs='?', const='7', default=None,
                         help='Select Centos Linux base image.')

parser.add_argument('--cuda', type=str, nargs='?', const='11.0', default=None,
                    help='Select a CUDA version for a base Linux image from NVIDIA.')

parser.add_argument('--mpi', type=str, nargs='?', const='openmpi', default=None,
                    help='Enable MPI (default disabled) and optionally select distribution (default: openmpi)')

parser.add_argument('--tsan', type=str, nargs='?', const='llvm', default=None,
                    help='Build special compiler versions with TSAN OpenMP support')

parser.add_argument('--hipsycl', type=str, nargs='?', default=None,
                    help='Select hipSYCL repository tag/commit/branch.')

parser.add_argument('--rocm', type=str, nargs='?', const='debian', default=None,
                    help='Select AMD compute engine version.')

parser.add_argument('--intel-compute-runtime', action='store_true', default=False,
                    help='Include Intel Compute Runtime.')

parser.add_argument('--clfft', type=str, nargs='?', const='master', default=None,
                    help='Add external clFFT libraries to the build image')

parser.add_argument('--heffte', type=str, nargs='?', default=None,
                    help='Select heffte repository tag/commit/branch.')

parser.add_argument('--doxygen', type=str, nargs='?', const='1.8.5', default=None,
                    help='Add doxygen environment for documentation builds. Also adds other requirements needed for final docs images.')

# Supported Python versions for maintained branches.
_python_versions = ['3.7.7', '3.8.2', '3.9.1']
parser.add_argument('--venvs', nargs='*', type=str, default=_python_versions,
                    help='List of Python versions ("major.minor.patch") for which to install venvs. '
                         'Default: {}'.format(' '.join(_python_versions)))


def image_name(configuration: argparse.Namespace) -> str:
    """Generate docker image name.

    Image names have the form ``ci-<slug>``, where the configuration slug has the form::

        <distro>-<version>-<compiler>-<major version>[-<gpusdk>-<version>][-<use case>]

    This function also applies an appropriate Docker image repository prefix.

    Arguments:
        configuration: Docker image configuration as described by the parsed arguments.

    """
    elements = []
    for distro in ('centos', 'ubuntu'):
        version = getattr(configuration, distro, None)
        if version is not None:
            elements.append(distro + '-' + version)
            break
    for compiler in ('llvm', 'gcc'):
        version = getattr(configuration, compiler, None)
        if version is not None:
            elements.append(compiler + '-' + str(version).split('.')[0])
            break
    for gpusdk in ('cuda', 'hipsycl'):
        version = getattr(configuration, gpusdk, None)
        if version is not None:
            elements.append(gpusdk + '-' + version)
    if configuration.oneapi is not None:
        elements.append('oneapi-' + configuration.oneapi)
    if configuration.intel_compute_runtime:
        elements.append('intel-compute-runtime')
    if configuration.rocm is not None:
        if (configuration.rocm != 'debian'):
            elements.append('rocm-' + configuration.rocm)

    # Check for special cases
    # The following attribute keys indicate the image is built for the named
    # special use case.
    cases = {'doxygen': 'docs',
             'tsan': 'tsan'}
    for attr in cases:
        value = getattr(configuration, attr, None)
        if value is not None:
            elements.append(cases[attr])
    slug = '-'.join(elements)
    # we are using the GitLab container registry to store the images
    # to get around issues with pulling them repeatedly from DockerHub
    # and running into the image pull limitation there.
    return 'registry.gitlab.com/gromacs/gromacs/ci-' + slug


if __name__ == "__main__":
    args = argparse.ArgumentParser(parents=[parser]).parse_args()
    print(image_name(args))
