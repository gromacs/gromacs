#!/usr/bin/env python
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

"""Building block based Dockerfile generation for CI testing images.

Generates a set of docker images used for running GROMACS CI on Gitlab.
The images are prepared according to a selection of build configuration targets
that hope to cover a broad enough scope of different possible systems,
allowing us to check compiler types and versions, as well as libraries used
for accelerators and parallel communication systems. Each combinations is
described as an entry in the build_configs dictionary, with the script
analysing the logic and adding build stages as needed.

Based on the example script provided by the NVidia HPCCM repository.

Reference:
    `NVidia HPC Container Maker <https://github.com/NVIDIA/hpc-container-maker>`__

Authors:
    * Paul Bauer <paul.bauer.q@gmail.com>
    * Eric Irrgang <ericirrgang@gmail.com>
    * Joe Jordan <e.jjordan12@gmail.com>
    * Mark Abraham <mark.j.abraham@gmail.com>
    * Gaurav Garg <gaugarg@nvidia.com>

Usage::

    $ python3 scripted_gmx_docker_builds.py --help
    $ python3 scripted_gmx_docker_builds.py --format docker > Dockerfile && docker build .
    $ python3 scripted_gmx_docker_builds.py | docker build -

See Also:
    :file:`buildall.sh`

"""

import argparse
import collections
import collections.abc
import copy
import os
import packaging.version
import shlex
import typing

try:
    import utility
except ImportError:
    raise RuntimeError(
        "This module assumes availability of supporting modules in the same directory. Add the directory to "
        "PYTHONPATH or invoke Python from within the module directory so module location can be resolved."
    )


def shlex_join(split_command):
    """Return a shell-escaped string from *split_command*.

    Copied from Python 3.8.
    Can be replaced with shlex.join once we don't need to support Python 3.7.
    """
    return " ".join(shlex.quote(arg) for arg in split_command)


# Basic packages for all final images.
_common_packages = [
    "build-essential",
    "ca-certificates",
    "ccache",
    "cmake",
    "git",
    "gnupg",
    "gpg-agent",
    "less",
    "libfftw3-dev",
    "libhwloc-dev",
    "liblapack-dev",
    "libx11-dev",
    "moreutils",
    "ninja-build",
    "rsync",
    "valgrind",
    "vim",
    "wget",
    "xsltproc",
]

_opencl_extra_packages = [
    "nvidia-opencl-dev",
    # The following require apt_ppas=['ppa:intel-opencl/intel-opencl'] on Ubuntu prior to 22.04
    "intel-opencl-icd",
    "ocl-icd-libopencl1",
    "ocl-icd-opencl-dev",
    "opencl-headers",
]

_rocm_extra_packages = [
    # The following require
    #             apt_keys=['http://repo.radeon.com/rocm/rocm.gpg.key'],
    #             apt_repositories=['deb [arch=amd64] http://repo.radeon.com/rocm/apt/X.Y.Z/ ubuntu main']
    "clinfo",
    "libelf1",
]

_rocm_version_dependent_packages = [
    # The following require
    #             apt_keys=['http://repo.radeon.com/rocm/rocm.gpg.key'],
    #             apt_repositories=['deb [arch=amd64] http://repo.radeon.com/rocm/apt/X.Y.Z/ ubuntu main']
    "hipfft",
    "hipfft-dev",
    "rocfft",
    "rocfft-dev",
    "rocprim",
    "rocprim-dev",
    "rocm-dev",
    "rocm-opencl",
]

# Extra packages required to build CP2K
_cp2k_extra_packages = [
    "autoconf",
    "autogen",
    "automake",
    "autotools-dev",
    "bzip2",
    "less",
    "libtool",
    "make",
    "nano",
    "patch",
    "pkg-config",
    "python",
    "python-numpy",
    "python3",
    "unzip",
    "xxd",
    "zlib1g-dev",
    "libopenblas-dev",
]

# Extra packages needed to build Intel Compute Runtime
_intel_compute_runtime_extra_packages = [
    "intel-opencl-icd",
    "intel-level-zero-gpu",
    "level-zero",
    "libmfx1",
]

# Extra packages needed to build Python installations from source.
# For Ubuntu 22.04, 'python-*' is later replaced with 'python3-*'
_python_extra_packages = [
    "build-essential",
    "ca-certificates",
    "ccache",
    "curl",
    "git",
    "libbz2-dev",
    "libffi-dev",
    "liblzma-dev",
    "libncurses5-dev",
    "libncursesw5-dev",
    "libreadline-dev",
    "libsqlite3-dev",
    "libssl-dev",
    "llvm",
    "python-openssl",
    "vim",
    "wget",
    "zlib1g-dev",
]

# Extra packages needed for images for building documentation.
_docs_extra_packages = [
    "autoconf",
    "automake",
    "autopoint",
    "autotools-dev",
    "bison",
    "flex",
    "ghostscript",
    "graphviz",
    "help2man",
    "imagemagick",
    "libtool",
    "mscgen",
    "m4",
    "openssh-client",
    "plantuml",
    "texinfo",
    "texlive-latex-base",
    "texlive-latex-extra",
    "texlive-fonts-recommended",
    "texlive-fonts-extra",
    "tex-gyre",
]

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="GROMACS CI image creation script", parents=[utility.parser]
)

parser.add_argument(
    "--format",
    type=str,
    default="docker",
    choices=["docker", "singularity"],
    help="Container specification format (default: %(default)s)",
)


def base_image_tag(args) -> str:
    """Generate *image* for hpccm.baseimage()."""
    # Check if we use CUDA images or plain linux images
    if args.cuda is not None:
        cuda_version_tag = "nvidia/cuda:" + args.cuda + "-devel"
        if args.centos is not None:
            cuda_version_tag += "-centos" + args.centos
        elif args.ubuntu is not None:
            cuda_version_tag += "-ubuntu" + args.ubuntu
        else:
            raise RuntimeError("Logic error: no Linux distribution selected.")

        base_image_tag = cuda_version_tag
    else:
        if args.centos is not None:
            base_image_tag = "centos:centos" + args.centos
        elif args.ubuntu is not None:
            base_image_tag = "ubuntu:" + args.ubuntu
        else:
            raise RuntimeError("Logic error: no Linux distribution selected.")
    return base_image_tag


def hpccm_distro_name(args) -> str:
    """Generate *_distro* for hpccm.baseimage().

    Convert the linux distribution variables into something that hpccm
    understands.

    The same format is used by the lower level hpccm.config.set_linux_distro().
    """
    if args.centos is not None:
        name_mapping = {"7": "centos7", "8": "centos8"}
        if args.centos in name_mapping:
            hpccm_name = name_mapping[args.centos]
        else:
            raise RuntimeError("Logic error: unsupported CentOS distribution selected.")
    elif args.ubuntu is not None:
        name_mapping = {
            "24.04": "ubuntu24",
            "22.04": "ubuntu22",
            "20.04": "ubuntu20",
            "18.04": "ubuntu18",
            "16.04": "ubuntu16",
        }
        if args.ubuntu in name_mapping:
            hpccm_name = name_mapping[args.ubuntu]
        else:
            raise RuntimeError("Logic error: unsupported Ubuntu distribution selected.")
    else:
        raise RuntimeError("Logic error: no Linux distribution selected.")
    return hpccm_name


def get_llvm_packages(args) -> typing.Iterable[str]:
    # If we use the package version of LLVM, we need to install extra packages for it.
    if (args.llvm is not None) and (args.tsan is None):
        packages = [
            f"libomp-{args.llvm}-dev",
            f"libomp5-{args.llvm}",
            "clang-format-" + str(args.llvm),
            "clang-tidy-" + str(args.llvm),
        ]
        if args.adaptivecpp is not None:
            packages += [
                f"llvm-{args.llvm}-dev",
                f"libclang-{args.llvm}-dev",
                f"libclang-rt-{args.llvm}-dev",
                f"lld-{args.llvm}",
            ]
        return packages
    else:
        return []


def get_opencl_packages(args) -> typing.List[str]:
    if (args.doxygen is None) and (args.oneapi is None):
        return _opencl_extra_packages
    else:
        return []


def get_rocm_packages(args) -> typing.List[str]:
    if args.rocm is None:
        return []
    else:
        packages = _rocm_extra_packages
        packages.extend(get_rocm_version_dependent_packages(args))
        return packages


def get_rocm_version_dependent_packages(args) -> typing.List[str]:
    packages = []
    rocm_version = args.rocm
    rocm_version_parsed = [int(i) for i in rocm_version.split(".")]
    # if the version does not contain the final patch version, add it manually to make sure that
    # we can install the packages
    if len(rocm_version_parsed) < 3:
        rocm_version = rocm_version + ".0"
    for entry in _rocm_version_dependent_packages:
        packages.append(entry + rocm_version)

    return packages


def get_rocm_repository(args) -> "hpccm.building_blocks.base":
    if args.ubuntu is None:
        raise RuntimeError("ROCm only supported on Ubuntu")
    else:
        rocm_version = [int(i) for i in args.rocm.split(".")]
        if rocm_version[0] < 5 or (rocm_version[0] == 5 and rocm_version[1] < 3):
            dist_string = "ubuntu"
        else:
            dist_string = {"20.04": "focal", "22.04": "jammy", "24.04": "noble"}[
                args.ubuntu
            ]
    return hpccm.building_blocks.packages(
        apt_keys=["http://repo.radeon.com/rocm/rocm.gpg.key"],
        apt_repositories=[
            f"deb [arch=amd64] http://repo.radeon.com/rocm/apt/{args.rocm}/ {dist_string} main"
        ],
    )


def get_cp2k_packages(args) -> typing.List[str]:
    cp2k_packages = []
    if args.cp2k:
        cp2k_packages.extend(_cp2k_extra_packages)
        if args.mpi is not None:
            cp2k_packages.append("libfftw3-mpi-dev")
    return cp2k_packages


def get_compiler(
    args, compiler_build_stage: "hpccm.Stage" = None
) -> "hpccm.building_blocks.base":
    # Compiler
    if args.llvm is not None:
        # Build our own version instead to get TSAN + OMP
        if args.tsan is not None:
            if compiler_build_stage is not None:
                compiler = compiler_build_stage.runtime(_from="tsan")
            else:
                raise RuntimeError("No TSAN compiler build stage!")
        # Build the default compiler if we don't need special support
        else:
            # Always use the "upstream" llvm repositories because the
            # main ubuntu repositories stop adding support for new
            # llvm versions after a few llvm releases.
            compiler = hpccm.building_blocks.llvm(version=args.llvm, upstream=True)

    elif args.oneapi is not None:
        if compiler_build_stage is not None:
            compiler = compiler_build_stage.runtime(_from="oneapi")
            # Prepare the toolchain (needed only for builds done within the Dockerfile, e.g.
            # OpenMPI builds, which don't currently work for other reasons)
            if args.oneapi.startswith("2024.0"):
                path = f"/opt/intel/oneapi/compiler/{args.oneapi}/linux/bin"
            else:
                path = f"/opt/intel/oneapi/compiler/{args.oneapi}/bin"
            oneapi_toolchain = hpccm.toolchain(
                CC=f"{path}/icx",
                CXX=f"{path}/icpx",
            )
            setattr(compiler, "toolchain", oneapi_toolchain)

        else:
            raise RuntimeError("No oneAPI compiler build stage!")

    elif args.intel_llvm is not None:
        if compiler_build_stage is not None:
            compiler = compiler_build_stage.runtime(_from="intel_llvm")
            intel_llvm_toolchain = hpccm.toolchain(
                CC="/opt/intel-llvm/bin/clang", CXX="/opt/intel-llvm/bin/clang++"
            )
            setattr(compiler, "toolchain", intel_llvm_toolchain)
        else:
            raise RuntimeError("No IntelLLVM compiler build stage!")

    elif args.gcc is not None:
        if args.cp2k is not None:
            compiler = hpccm.building_blocks.gnu(
                extra_repository=True, version=args.gcc, fortran=True
            )
        else:
            compiler = hpccm.building_blocks.gnu(
                extra_repository=True, version=args.gcc, fortran=False
            )
    else:
        raise RuntimeError("Logic error: no compiler toolchain selected.")
    return compiler


def get_gdrcopy(args, compiler):
    if args.cuda is not None:
        if hasattr(compiler, "toolchain"):
            # Version last updated August 15, 2024
            return hpccm.building_blocks.gdrcopy(
                toolchain=compiler.toolchain, version="2.4.1"
            )
        else:
            raise RuntimeError("compiler is not an HPCCM compiler building block!")
    else:
        return None


def get_ucx(args, compiler, gdrcopy):
    if args.cuda is not None or args.rocm is not None:
        if hasattr(compiler, "toolchain"):
            use_gdrcopy = gdrcopy is not None
            # We disable `-Werror`, since there are some unknown pragmas and unused variables which upset clang
            toolchain = copy.copy(compiler.toolchain)
            toolchain.CFLAGS = "-Wno-error"
            configure_opts = []
            if args.rocm is not None:
                configure_opts.append("--with-rocm=/opt/rocm")
            if args.cuda is not None:
                configure_opts.append("--with-cuda=/usr/local/cuda")
            # Version last updated August 15, 2024
            return hpccm.building_blocks.ucx(
                toolchain=toolchain,
                gdrcopy=use_gdrcopy,
                version="1.17.0",
                configure_opts=configure_opts,
            )
        else:
            raise RuntimeError("compiler is not an HPCCM compiler building block!")
    else:
        return None


def get_mpi(args, compiler, ucx):
    # If needed, add MPI to the image
    if args.mpi is not None:
        if args.mpi == "openmpi":
            if hasattr(compiler, "toolchain"):
                if args.oneapi is not None:
                    raise RuntimeError("oneAPI building OpenMPI is not supported")
                use_cuda = args.cuda is not None
                use_ucx = ucx is not None
                # Version last updated October 7, 2022
                return hpccm.building_blocks.openmpi(
                    toolchain=compiler.toolchain,
                    version="4.1.4",
                    cuda=use_cuda,
                    ucx=use_ucx,
                    infiniband=False,
                )
            else:
                raise RuntimeError("compiler is not an HPCCM compiler building block!")
        elif args.mpi == "mpich":
            if hasattr(compiler, "toolchain"):
                use_cuda = args.cuda is not None
                use_rocm = args.rocm is not None
                use_ucx = ucx is not None
                flags = {}
                if ucx is not None:
                    flags["with-device"] = "ch4:ucx"
                mpich_stage = hpccm.Stage()
                # Python needed for configuring
                mpich_stage += hpccm.building_blocks.python(
                    python3=True, python2=False, devel=False
                )
                # Version last updated August 15, 2024
                mpich_stage += hpccm.building_blocks.mpich(
                    toolchain=compiler.toolchain,
                    version="4.2.2",
                    cuda=use_cuda,
                    rocm=use_rocm,
                    ucx=use_ucx,
                    infiniband=False,
                    disable_fortran=True,
                    **flags,
                )
                return mpich_stage
            else:
                raise RuntimeError("compiler is not an HPCCM compiler building block!")
        elif args.mpi == "impi":
            # TODO Intel MPI from the oneAPI repo is not working reliably,
            # reasons are unclear. When solved, add packagages called:
            # 'intel-oneapi-mpi', 'intel-oneapi-mpi-devel'
            # during the compiler stage.
            # TODO also consider hpccm's intel_mpi package if that doesn't need
            # a license to run.
            raise RuntimeError("Intel MPI recipe not implemented yet.")
        else:
            raise RuntimeError("Requested unknown MPI implementation.")
    else:
        return None


def get_oneapi_plugins(args):
    # To get this token, register at https://developer.codeplay.com/ and generate new API token on the "Setting" page.
    # Then place the toke in this environment variable when building the container.
    token = os.getenv("CODEPLAY_API_TOKEN")
    blocks = []

    def _add_plugin(variant):
        if args.oneapi is None:
            raise RuntimeError("Cannot install oneAPI plugins without oneAPI.")
        if token is None:
            raise RuntimeError(
                "Need CODEPLAY_API_TOKEN env. variable to install oneAPI plugins"
            )
        backend_version = {"nvidia": args.cuda, "amd": args.rocm}[variant]
        if backend_version.count(".") == 2:
            backend_version = ".".join(backend_version.split(".")[:2])  # 12.0.1 -> 12.0
        oneapi_version = args.oneapi
        url = f"https://developer.codeplay.com/api/v1/products/download?product=oneapi&version={oneapi_version}&variant={variant}&filters[]=linux&filters[]={backend_version}&aat={token}"
        outfile = f"/tmp/oneapi_plugin_{variant}.sh"
        blocks.append(
            hpccm.primitives.shell(
                commands=[
                    f"wget --content-disposition '{url}' --output-document '{outfile}'",
                    f"bash '{outfile}' --yes",
                ]
            )
        )

    if args.oneapi_plugin_nvidia:
        _add_plugin("nvidia")
    if args.oneapi_plugin_amd:
        _add_plugin("amd")
        # Need to make ROCm libraries discoverable
        blocks.append(
            hpccm.primitives.shell(
                commands=[
                    "echo '/opt/rocm/lib/' > /etc/ld.so.conf.d/rocm.conf",
                    "ldconfig",
                ]
            )
        )
    return blocks


def get_clfft(args):
    if args.clfft is not None:
        return hpccm.building_blocks.generic_cmake(
            repository="https://github.com/clMathLibraries/clFFT.git",
            prefix="/usr/local",
            recursive=True,
            branch=args.clfft,
            directory="clFFT/src",
        )
    else:
        return None


def get_heffte(args):
    if args.heffte is not None:
        return hpccm.building_blocks.generic_cmake(
            cmake_opts=[
                "-D CMAKE_BUILD_TYPE=Release",
                "-D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda",
                "-D Heffte_ENABLE_CUDA=ON",
                "-D Heffte_ENABLE_FFTW=OFF",
                "-D BUILD_SHARED_LIBS=ON",
            ],
            repository="https://bitbucket.org/icl/heffte.git",
            prefix="/usr/local",
            recursive=True,
            commit=args.heffte,
            directory="heffte",
        )
    else:
        return None


def get_nvhpcsdk(args):
    if args.nvhpcsdk is not None:
        return hpccm.building_blocks.nvhpc(
            eula=True,
            cuda_multi=False,
            environment=False,
            mpi=False,
            version=args.nvhpcsdk,
        )
    else:
        return None


def get_adaptivecpp(args):
    if args.adaptivecpp is None:
        return None
    if args.rocm is None:
        raise RuntimeError("AdaptiveCpp requires the ROCm packages")
    if args.llvm is None:
        # We're using ROCm LLVM in this case, which is not compatible with CUDA
        if args.cuda is not None:
            raise RuntimeError(
                "Can not build AdaptiveCpp with CUDA and no upstream LLVM"
            )
    preconfigure = []

    if args.llvm is not None:
        cmake_opts = [
            "-DCMAKE_C_COMPILER=clang-{}".format(args.llvm),
            "-DCMAKE_CXX_COMPILER=clang++-{}".format(args.llvm),
            "-DLLVM_DIR=/usr/lib/llvm-{}/cmake/".format(args.llvm),
        ]
        if int(args.llvm) >= 18:  # and args.adaptivecpp <= 24.06:
            preconfigure.append(
                "ln -s /usr/lib/llvm-{0}/lib/libLLVM-{0}.so /usr/lib/llvm-{0}/lib/libLLVM.so".format(
                    args.llvm
                )
            )
    else:
        cmake_opts = [
            "-DCMAKE_C_COMPILER=/opt/rocm/bin/amdclang",
            "-DCMAKE_CXX_COMPILER=/opt/rocm/bin/amdclang++",
            "-DLLVM_DIR=/opt/rocm/llvm/lib/cmake/llvm",
            "-DWITH_SSCP_COMPILER=OFF",
            "-DWITH_OPENCL_BACKEND=OFF",
            "-DWITH_LEVEL_ZERO_BACKEND=OFF",
        ]

    cmake_opts += [
        "-DCMAKE_PREFIX_PATH=/opt/rocm/lib/cmake",
        "-DWITH_ROCM_BACKEND=ON",
    ]

    if args.cuda is not None:
        cmake_opts += [
            "-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda",
            "-DWITH_CUDA_BACKEND=ON",
        ]

    adaptivecpp_version_opts = {}
    if "." in args.adaptivecpp:
        adaptivecpp_version_opts["branch"] = "v" + args.adaptivecpp
    else:
        adaptivecpp_version_opts["commit"] = args.adaptivecpp

    return hpccm.building_blocks.generic_cmake(
        repository="https://github.com/AdaptiveCpp/AdaptiveCpp.git",
        directory="/var/tmp/AdaptiveCpp",
        prefix="/usr/local",
        recursive=True,
        preconfigure=preconfigure,
        cmake_opts=["-DCMAKE_BUILD_TYPE=Release", *cmake_opts],
        **adaptivecpp_version_opts,
    )


def get_cp2k(args):
    if args.cp2k is None:
        return None

    if args.gcc is None:
        raise RuntimeError("CP2K build requires GNU compilers")

    make_commands = ["make -j$(nproc) ARCH=local VERSION=ssmp libcp2k"]
    if args.mpi is not None:
        make_commands += ["make -j$(nproc) ARCH=local VERSION=psmp libcp2k"]
    make_commands += ["rm -rf ./obj"]

    return hpccm.building_blocks.generic_build(
        repository="https://github.com/cp2k/cp2k.git",
        branch=f"support/v{args.cp2k}",
        recursive=True,
        build=[
            "mkdir -p /opt/cp2k",
            "cp -rf ./* /opt/cp2k/",
            "cd /opt/cp2k/tools/toolchain",
            "./install_cp2k_toolchain.sh  \
                --with-gcc=system            \
                --with-cmake=system          \
                --with-openmpi=system        \
                --with-fftw=system           \
                --with-openblas=system       \
                --with-scalapack=install     \
                --with-gsl=no                \
                --with-elpa=no               \
                --with-spglib=no             \
                --with-spfft=no              \
                --with-cosma=no              \
                --with-libvori=no            \
                --with-sirius=no             \
                --with-hdf5=no               \
                --with-libxc=install         \
                --with-libxsmm=install       \
                --with-libint=install        \
                --libint-lmax=5",
            "rm -rf ./build",
            "cp ./install/arch/local.* ../../arch/",
            "cd ../../",
        ]
        + make_commands,
    )


def add_tsan_compiler_build_stage(
    input_args, output_stages: typing.Mapping[str, "hpccm.Stage"]
):
    """Isolate the expensive TSAN preparation stage.

    This is a very expensive stage, but has few and disjoint dependencies, and
    its output is easily compartmentalized (/usr/local) so we can isolate this
    build stage to maximize build cache hits and reduce rebuild time, bookkeeping,
    and final image size.
    """
    if not isinstance(output_stages, collections.abc.MutableMapping):
        raise RuntimeError("Need output_stages container.")
    if "compiler_build" in output_stages:
        raise RuntimeError('"compiler_build" output stage is already present.')
    tsan_stage = hpccm.Stage()
    tsan_stage += hpccm.primitives.baseimage(
        image=base_image_tag(input_args),
        _distro=hpccm_distro_name(input_args),
        _as="tsan",
    )

    tsan_stage += hpccm.building_blocks.packages(
        ospackages=["git", "ca-certificates", "build-essential", "cmake"]
    )
    # CMake will get duplicated later, but this is an expensive image, and it isn't worth optimizing
    # out that duplication...
    tsan_stage += hpccm.building_blocks.python(python3=True, python2=False, devel=False)

    compiler_branch = "release/" + str(input_args.llvm) + ".x"
    tsan_stage += hpccm.building_blocks.generic_cmake(
        repository="https://github.com/llvm/llvm-project.git",
        directory="/var/tmp/llvm-project/llvm/",
        prefix="/usr/local",
        recursive=True,
        branch=compiler_branch,
        cmake_opts=[
            "-D CMAKE_BUILD_TYPE=Release",
            '-D LLVM_ENABLE_PROJECTS="clang;openmp;clang-tools-extra;compiler-rt;lld"',
            "-D LIBOMP_TSAN_SUPPORT=on",
        ],
        postinstall=[
            "ln -s /usr/local/bin/clang++ /usr/local/bin/clang++-"
            + str(input_args.llvm),
            "ln -s /usr/local/bin/clang-format /usr/local/bin/clang-format-"
            + str(input_args.llvm),
            "ln -s /usr/local/bin/clang-tidy /usr/local/bin/clang-tidy-"
            + str(input_args.llvm),
            "ln -s /usr/local/share/clang/run-clang-tidy.py /usr/local/bin/run-clang-tidy-"
            + str(input_args.llvm)
            + ".py",
            "ln -s /usr/local/bin/run-clang-tidy-"
            + str(input_args.llvm)
            + ".py /usr/local/bin/run-clang-tidy-"
            + str(input_args.llvm),
            "ln -s /usr/local/libexec/c++-analyzer /usr/local/bin/c++-analyzer-"
            + str(input_args.llvm),
        ],
    )
    output_stages["compiler_build"] = tsan_stage


def oneapi_runtime(_from="0"):
    oneapi_runtime_stage = hpccm.Stage()
    oneapi_runtime_stage += hpccm.primitives.copy(
        _from="oneapi-build",
        files={"/opt/intel": "/opt/intel", "/etc/bash.bashrc": "/etc/bash.bashrc"},
    )
    return oneapi_runtime_stage


def add_oneapi_compiler_build_stage(
    input_args, output_stages: typing.Mapping[str, "hpccm.Stage"]
):
    """Isolate the oneAPI preparation stage.

    This stage is isolated so that its installed components are minimized in the
    final image (chiefly /opt/intel) and its environment setup script can be
    sourced. This also helps with rebuild time and final image size.
    """
    if not isinstance(output_stages, collections.abc.MutableMapping):
        raise RuntimeError("Need output_stages container.")
    if "compiler_build" in output_stages:
        raise RuntimeError('"compiler_build" output stage is already present.')
    oneapi_stage = hpccm.Stage()
    oneapi_stage += hpccm.primitives.baseimage(
        image=base_image_tag(input_args),
        _distro=hpccm_distro_name(input_args),
        _as="oneapi-build",
    )

    version = str(input_args.oneapi)

    # Add required components for the next stage (both for hpccm and Intel's setvars.sh script)
    oneapi_stage += hpccm.building_blocks.packages(
        ospackages=["wget", "gnupg2", "ca-certificates", "lsb-release"]
    )
    oneapi_stage += hpccm.building_blocks.packages(
        apt_keys=[
            "https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB"
        ],
        apt_repositories=["deb https://apt.repos.intel.com/oneapi all main"],
        # Add minimal packages (not the whole HPC toolkit!)
        ospackages=[
            f"intel-oneapi-dpcpp-cpp-{version}",
            f"intel-oneapi-openmp-{version}",
            f"intel-oneapi-mkl-{version}",
            f"intel-oneapi-mkl-devel-{version}",
        ],
    )
    # Ensure that all bash shells on the final container will have access to oneAPI
    oneapi_stage += hpccm.primitives.shell(
        commands=[
            'echo "source /opt/intel/oneapi/setvars.sh" >> /etc/bash.bashrc',
            "unlink /opt/intel/oneapi/compiler/latest",
            f"ln -sf /opt/intel/oneapi/compiler/{version} /opt/intel/oneapi/compiler/latest",
        ]
    )
    setattr(oneapi_stage, "runtime", oneapi_runtime)

    output_stages["compiler_build"] = oneapi_stage


def intel_llvm_runtime(_from="0"):
    llvm_runtime_stage = hpccm.Stage()
    llvm_runtime_stage += hpccm.primitives.copy(
        _from="intel-llvm-build", files={"/opt/intel-llvm": "/opt/intel-llvm"}
    )

    bashrc = [
        "export DPCPP_HOME=/opt/intel-llvm",
        "export PATH=${DPCPP_HOME}/bin:$PATH",
        "export LD_LIBRARY_PATH=${DPCPP_HOME}/lib:${LD_LIBRARY_PATH}",
        'export CFLAGS="-isystem ${DPCPP_HOME}/include"',
        'export CXXFLAGS="-isystem ${DPCPP_HOME}/include"',
    ]
    # Since we cannot just create a file, we write to it line-by-line using "echo".
    # We must shlex.quote the lines to ensure all spaces/quotes/etc are preserved.
    commands = [
        "echo {} >> /opt/intel-llvm/setenv.sh".format(shlex.quote(line))
        for line in bashrc
    ] + ['echo source "/opt/intel-llvm/setenv.sh" >> /etc/bash.bashrc']
    llvm_runtime_stage += hpccm.primitives.shell(commands=commands)

    return llvm_runtime_stage


def add_intel_llvm_compiler_build_stage(
    input_args, output_stages: typing.Mapping[str, "hpccm.Stage"]
):
    """Isolate the Intel LLVM (open-source oneAPI) preparation stage.

    This stage is isolated so that its installed components are minimized in the
    final image (chiefly /opt/intel) and its environment setup script can be
    sourced. This also helps with rebuild time and final image size.
    """
    if not isinstance(output_stages, collections.abc.MutableMapping):
        raise RuntimeError("Need output_stages container.")
    if "compiler_build" in output_stages:
        raise RuntimeError('"compiler_build" output stage is already present.')
    llvm_stage = hpccm.Stage()
    llvm_stage += hpccm.primitives.baseimage(
        image=base_image_tag(input_args),
        _distro=hpccm_distro_name(input_args),
        _as="intel-llvm-build",
    )
    buildbot_flags = [
        "--build-type=Release",
        "--cuda",  # Build with CUDA support
        "--llvm-external-projects=openmp",  # Enable OpenMP
        "--obj-dir=/var/tmp/llvm/llvm/build",  # Build directory
        # Disable OpenMP offload targets we're not using
        "--cmake-opt=-DLIBOMPTARGET_BUILD_AMDGPU_PLUGIN=FALSE",
        "--cmake-opt=-DLIBOMPTARGET_BUILD_CUDA_PLUGIN=FALSE",
        "--cmake-opt=-DLIBOMPTARGET_BUILD_DEVICERTL_BCLIB=FALSE",
        # Help CMake find CUDA Driver stub, see https://github.com/opencv/opencv/issues/6577
        "--cmake-opt=-DCMAKE_LIBRARY_PATH=/usr/local/cuda/targets/x86_64-linux/lib/stubs/",
    ]

    llvm_stage += hpccm.building_blocks.packages(
        ospackages=[
            "git",
            "ninja-build",
            "cmake",
            "python3",
            "python3-dev",
            "build-essential",
            "wget",
        ]
    )

    if args.rocm is not None:
        llvm_stage += get_rocm_repository(args)
        llvm_stage += hpccm.building_blocks.packages(
            ospackages=get_rocm_packages(args), aptitude=True
        )
        buildbot_flags.extend(["--hip", "--hip-platform", "AMD"])

    llvm_stage += hpccm.building_blocks.generic_build(
        repository="https://github.com/intel/llvm.git",
        directory="llvm/llvm",
        build=[
            "mkdir -p /var/tmp/llvm/llvm/build",
            shlex_join(
                ["python3", "/var/tmp/llvm/buildbot/configure.py", *buildbot_flags]
            ),
            "cd /var/tmp/llvm/llvm/build",
            # Must be called after the configure.py
            shlex_join(
                [
                    "cmake",
                    "/var/tmp/llvm/llvm",
                    "-DCMAKE_INSTALL_PREFIX=/opt/intel-llvm/",
                ]
            ),
            "ninja",
            "ninja sycl-toolchain install",
        ],
        install=[],
        branch=input_args.intel_llvm,
    )

    setattr(llvm_stage, "runtime", intel_llvm_runtime)

    output_stages["compiler_build"] = llvm_stage


def prepare_venv(version: packaging.version.Version) -> typing.Sequence[str]:
    """Get shell commands to set up the venv for the requested Python version."""
    major = version.major
    minor = version.minor  # type: int

    pyenv = "$HOME/.pyenv/bin/pyenv"

    py_ver = f"{major}.{minor}"
    venv_path = f"$HOME/venv/py{py_ver}"
    commands = [
        f"$({pyenv} prefix `{pyenv} whence python{py_ver}`)/bin/python -m venv {venv_path}"
    ]

    commands.append(f"{venv_path}/bin/python -m pip install --upgrade pip setuptools")
    # Install dependencies for building and testing gmxapi Python package.
    # WARNING: Please keep this list synchronized with python_packaging/gmxapi/requirements.txt
    #  and docs/requirements.txt
    # TODO: Get requirements.txt from an input argument.
    commands.append(
        f"""{venv_path}/bin/python -m pip install --upgrade \
            'black' \
            'breathe' \
            'build' \
            'cmake>=3.18.4' \
            'flake8>=3.7.7' \
            'furo' \
            'gcovr>=4.2' \
            'importlib-resources;python_version<"3.10"' \
            'mpi4py>=3.0.3' \
            'mypy' \
            'networkx>=2.0' \
            'numpy>1.7' \
            'packaging' \
            'pip>=10.1' \
            'pybind11>2.6' \
            'Pygments>=2.2.0' \
            'pytest>=4.6' \
            'python-gitlab' \
            'setuptools>=61' \
            'Sphinx>=4.0' \
            'sphinx-argparse' \
            'sphinx-copybutton' \
            'sphinx_inline_tabs' \
            'sphinxcontrib-autoprogram' \
            'sphinxcontrib-plantuml>=0.14' \
            'versioningit>=2' \
            'wheel'"""
    )
    return commands


def get_cmake_stages(*, input_args: argparse.Namespace, base: str):
    """Get the stage(s) necessary for the requested CMake versions.

    One (intermediate) build stage is created
    for each CMake version, based on the *base* stage.
    See ``--cmake`` option.

    Each stage uses the version number to determine an installation location:
        /usr/local/cmake-{version}

    The resulting path is easily copied into the main stage.

    Returns:
        dict of isolated CMake installation stages with keys from ``cmake-{version}``
    """
    cmake_stages = {}
    for cmake_version in input_args.cmake:
        stage_name = f"cmake-{cmake_version}"
        cmake_stages[stage_name] = hpccm.Stage()
        cmake_stages[stage_name] += hpccm.primitives.baseimage(
            image=base, _distro=hpccm_distro_name(input_args), _as=stage_name
        )
        cmake_stages[stage_name] += hpccm.building_blocks.cmake(
            eula=True, prefix=f"/usr/local/{stage_name}", version=cmake_version
        )
    return cmake_stages


def add_python_stages(
    input_args: argparse.Namespace,
    *,
    base: str,
    output_stages: typing.MutableMapping[str, "hpccm.Stage"],
):
    """Add the stage(s) necessary for the requested venvs.

    One intermediate build stage is created for each venv (see --venv option).

    Each stage partially populates Python installations and venvs in the home
    directory. The home directory is collected by the 'pyenv' stage for use by
    the main build stage.
    """
    if len(input_args.venvs) < 1:
        raise RuntimeError("No venvs to build...")
    if output_stages is None or not isinstance(output_stages, collections.abc.Mapping):
        raise RuntimeError("Need a container for output stages.")

    # Main Python stage that collects the environments from individual stages.
    # We collect the stages individually, rather than chaining them, because the
    # copy is a bit slow and wastes local Docker image space for each filesystem
    # layer.
    pyenv_stage = hpccm.Stage()
    pyenv_stage += hpccm.primitives.baseimage(
        image=base, _distro=hpccm_distro_name(input_args), _as="pyenv"
    )
    python_extra_packages = _python_extra_packages
    if input_args.ubuntu is not None and input_args.ubuntu != "20.04":
        python_extra_packages = [
            i.replace("python-", "python3-") for i in python_extra_packages
        ]
    pyenv_stage += hpccm.building_blocks.packages(ospackages=python_extra_packages)

    for version in [
        packaging.version.parse(py_ver) for py_ver in sorted(input_args.venvs)
    ]:
        stage_name = "py" + str(version)
        stage = hpccm.Stage()
        stage += hpccm.primitives.baseimage(
            image=base, _distro=hpccm_distro_name(input_args), _as=stage_name
        )
        stage += hpccm.building_blocks.packages(ospackages=python_extra_packages)

        # TODO: Use a non-root user for testing and Python virtual environments.
        stage += hpccm.primitives.shell(
            commands=[
                "curl https://pyenv.run | bash",
                """echo 'export PYENV_ROOT="$HOME/.pyenv"' >> $HOME/.bashrc""",
                """echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> $HOME/.bashrc""",
                """echo 'eval "$(pyenv init -)"' >> $HOME/.bashrc""",
                """echo 'eval "$(pyenv virtualenv-init -)"' >> $HOME/.bashrc""",
            ]
        )
        pyenv = "$HOME/.pyenv/bin/pyenv"
        commands = [
            f'PYTHON_CONFIGURE_OPTS="--enable-shared" {pyenv} install -s {version}'
        ]
        stage += hpccm.primitives.shell(commands=commands)

        commands = prepare_venv(version)
        stage += hpccm.primitives.shell(commands=commands)

        # TODO: Update user home directory.
        pyenv_stage += hpccm.primitives.copy(
            _from=stage_name, _mkdir=True, src=["/root/"], dest="/root"
        )

        # Add the intermediate build stage to the sequence
        output_stages[stage_name] = stage

    # TODO: If we activate pyenv for login shells, the `global` "version" should be full-featured.
    # # `version` should be a system installation or pyenv environment (or pyenv-virtualenv)
    # # with the dependencies for all of the Python aspects of CMake-driven builds.
    # commands = '{pyenv} global {version}'.format(
    #             pyenv=pyenv,
    #             version=...)
    # pyenv_stage += hpccm.primitives.shell(commands=commands)

    # Add the aggregating build stage to the sequence. This allows the main stage to copy
    # the files in a single stage, potentially reducing the overall output image size.
    output_stages["pyenv"] = pyenv_stage


def add_documentation_dependencies(
    input_args, output_stages: typing.MutableMapping[str, "hpccm.Stage"]
):
    """Add appropriate layers according to doxygen input arguments."""
    if input_args.doxygen is None:
        return
    # Always clone the same version of linkchecker (latest release at June 1, 2021)
    output_stages["main"] += hpccm.building_blocks.pip(
        pip="pip3",
        packages=["git+https://github.com/linkchecker/linkchecker.git@v10.0.1"],
    )
    output_stages["main"] += hpccm.primitives.shell(
        commands=[
            'sed -i \'/"XPS"/d;/"PDF"/d;/"PS"/d;/"EPS"/d;/disable ghostscript format types/d\' /etc/ImageMagick-6/policy.xml'
        ]
    )
    if input_args.doxygen == "1.8.5":
        doxygen_commit = "ed4ed873ab0e7f15116e2052119a6729d4589f7a"
        output_stages["main"] += hpccm.building_blocks.generic_autotools(
            repository="https://github.com/westes/flex.git",
            commit="f7788a9a0ecccdc953ed12043ccb59ca25714018",
            prefix="/tmp/install-of-flex",
            configure_opts=["--disable-shared"],
            preconfigure=["./autogen.sh"],
        )
        output_stages["main"] += hpccm.building_blocks.generic_autotools(
            repository="https://github.com/doxygen/doxygen.git",
            commit=doxygen_commit,
            prefix="",
            configure_opts=["--flex /tmp/install-of-flex/bin/flex", "--static"],
        )
    else:
        version = input_args.doxygen
        archive_name = f"doxygen-{version}.linux.bin.tar.gz"
        archive_url = f"https://sourceforge.net/projects/doxygen/files/rel-{version}/{archive_name}"
        binary_path = f"doxygen-{version}/bin/doxygen"
        commands = [
            "mkdir doxygen && cd doxygen",
            f"wget {archive_url}",
            f"tar xf {archive_name} {binary_path}",
            f"cp {binary_path} /usr/local/bin/",
            "cd .. && rm -rf doxygen",
        ]
        output_stages["main"] += hpccm.primitives.shell(commands=commands)


def add_base_stage(
    name: str, input_args, output_stages: typing.MutableMapping[str, "hpccm.Stage"]
):
    """Establish dependencies that are shared by multiple parallel stages."""
    # Building blocks are chunks of container-builder instructions that can be
    # copied to any build stage with the addition operator.
    building_blocks = collections.OrderedDict()
    building_blocks["base_packages"] = hpccm.building_blocks.packages(
        ospackages=_common_packages
    )

    # These are the most expensive and most reusable layers, so we put them first.
    building_blocks["compiler"] = get_compiler(
        input_args, compiler_build_stage=output_stages.get("compiler_build")
    )
    if args.rocm is not None:
        building_blocks["rocm"] = [
            get_rocm_repository(args),
            hpccm.building_blocks.packages(ospackages=get_rocm_packages(args)),
        ]
    building_blocks["gdrcopy"] = get_gdrcopy(input_args, building_blocks["compiler"])
    building_blocks["ucx"] = get_ucx(
        input_args, building_blocks["compiler"], building_blocks["gdrcopy"]
    )
    building_blocks["mpi"] = get_mpi(
        input_args, building_blocks["compiler"], building_blocks["ucx"]
    )

    # Create the stage from which the targeted image will be tagged.
    output_stages[name] = hpccm.Stage()

    output_stages[name] += hpccm.primitives.baseimage(
        image=base_image_tag(input_args),
        _distro=hpccm_distro_name(input_args),
        _as=name,
    )
    for bb in building_blocks.values():
        if bb is not None:
            output_stages[name] += bb


def build_stages(args) -> typing.Iterable["hpccm.Stage"]:
    """Define and sequence the stages for the recipe corresponding to *args*."""

    # A Dockerfile or Singularity recipe can have multiple build stages.
    # The main build stage can copy files from previous stages, though only
    # the last stage is included in the tagged output image. This means that
    # large or expensive sets of build instructions can be isolated in
    # local/temporary images, but all of the stages need to be output by this
    # script, and need to occur in the correct order, so we create a sequence
    # object early in this function.
    stages = collections.OrderedDict()

    # If we need TSAN or oneAPI support the early build is more complex,
    # so that our compiler images don't have all the cruft needed to get those things
    # installed.
    if args.llvm is not None and args.tsan is not None:
        add_tsan_compiler_build_stage(input_args=args, output_stages=stages)
    if args.oneapi is not None:
        add_oneapi_compiler_build_stage(input_args=args, output_stages=stages)
    if args.intel_llvm is not None:
        add_intel_llvm_compiler_build_stage(input_args=args, output_stages=stages)

    add_base_stage(name="build_base", input_args=args, output_stages=stages)

    # Add Python environments to MPI images, only, so we don't have to worry
    # about whether to install mpi4py.
    if args.mpi is not None and len(args.venvs) > 0:
        add_python_stages(base="build_base", input_args=args, output_stages=stages)

    # Building blocks are chunks of container-builder instructions that can be
    # copied to any build stage with the addition operator.
    building_blocks = collections.OrderedDict()

    # Install additional packages early in the build to optimize Docker build layer cache.
    os_packages = (
        list(get_llvm_packages(args))
        + get_opencl_packages(args)
        + get_cp2k_packages(args)
    )
    if args.doxygen is not None:
        os_packages += _docs_extra_packages
    if args.oneapi is not None:
        os_packages += ["lsb-release"]
    if args.adaptivecpp is not None:
        os_packages += ["libboost-fiber-dev"]
    building_blocks["extra_packages"] = []
    if args.intel_compute_runtime:
        repo = {
            "24.04": "deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu noble arc",
            "22.04": "deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu jammy arc",
            "20.04": "deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu focal main",
        }
        building_blocks["extra_packages"] += hpccm.building_blocks.packages(
            apt_keys=["https://repositories.intel.com/graphics/intel-graphics.key"],
            apt_repositories=[repo[args.ubuntu]],
        )
        os_packages += _intel_compute_runtime_extra_packages

    if args.ubuntu is not None and args.ubuntu == "20.04":
        building_blocks["extra_packages"] += hpccm.building_blocks.packages(
            ospackages=os_packages, apt_ppas=["ppa:intel-opencl/intel-opencl"]
        )
    else:
        # For 22.04, everything is packaged. Therefore, we don't need Intel repo, but we need to
        # use aptitude to resolve the conflicts between Ubuntu ROCm packages and AMD ROCm packages.
        building_blocks["extra_packages"] += hpccm.building_blocks.packages(
            ospackages=os_packages, aptitude=True
        )

    building_blocks["CP2K"] = get_cp2k(args)

    if args.cuda is not None and args.llvm is not None:
        # Hack to tell clang what version of CUDA we're using
        # based on https://github.com/llvm/llvm-project/blob/1fdec59bffc11ae37eb51a1b9869f0696bfd5312/clang/lib/Driver/ToolChains/Cuda.cpp#L43
        cuda_version_split = args.cuda.split(".")
        # LLVM requires having the version in x.y.z format, while args.cuda be be either x.y or x.y.z
        cuda_version_str = "{}.{}.{}".format(
            cuda_version_split[0],
            cuda_version_split[1],
            cuda_version_split[2] if len(cuda_version_split) > 2 else 0,
        )
        building_blocks["cuda-clang-workaround"] = hpccm.primitives.shell(
            commands=[
                f'echo "CUDA Version {cuda_version_str}" > /usr/local/cuda/version.txt'
            ]
        )

    building_blocks["oneapi_plugins"] = get_oneapi_plugins(args)

    building_blocks["clfft"] = get_clfft(args)

    building_blocks["heffte"] = get_heffte(args)

    building_blocks["nvhpcsdk"] = get_nvhpcsdk(args)
    if building_blocks["nvhpcsdk"] is not None:
        nvshmem_lib_path = (
            "/opt/nvidia/hpc_sdk/Linux_x86_64/"
            + args.nvhpcsdk
            + "/comm_libs/nvshmem/lib/:$LD_LIBRARY_PATH"
        )
        building_blocks["nvhpcsdk_path"] = hpccm.primitives.environment(
            variables={"LD_LIBRARY_PATH": nvshmem_lib_path}
        )

    building_blocks["AdaptiveCpp"] = get_adaptivecpp(args)

    # Add Python environments to MPI images, only, so we don't have to worry
    # about whether to install mpi4py.
    if args.mpi is not None and len(args.venvs) > 0:
        add_python_stages(base="build_base", input_args=args, output_stages=stages)

    cmake_stages = get_cmake_stages(input_args=args, base="build_base")
    stages.update(cmake_stages)

    # Create the stage from which the targeted image will be tagged.
    stages["main"] = hpccm.Stage()

    stages["main"] += hpccm.primitives.baseimage(
        image="build_base", _distro=hpccm_distro_name(args), _as="main"
    )
    for bb in building_blocks.values():
        if bb is not None:
            stages["main"] += bb

    # We always add Python3 and Pip
    stages["main"] += hpccm.building_blocks.python(python3=True, python2=False)

    # Add documentation requirements (doxygen and sphinx + misc).
    if args.doxygen is not None:
        add_documentation_dependencies(args, stages)

    for stage_name in cmake_stages:
        stages["main"] += hpccm.primitives.copy(
            _from=stage_name,
            _mkdir=True,
            src=[f"/usr/local/{stage_name}/"],
            dest=f"/usr/local/{stage_name}",
        )

    if "pyenv" in stages and stages["pyenv"] is not None:
        stages["main"] += hpccm.primitives.copy(
            _from="pyenv", _mkdir=True, src=["/root/.pyenv/"], dest="/root/.pyenv"
        )
        stages["main"] += hpccm.primitives.copy(
            _from="pyenv", _mkdir=True, src=["/root/venv/"], dest="/root/venv"
        )
        # TODO: Update user home directory.
        # TODO: If we activate pyenv for login shells, the `global` "version" should be full-featured.
        # stages['main'] += hpccm.primitives.copy(_from='pyenv', src=['/root/.bashrc'],
        #                                         dest='/root/')

    # Make sure that `python` resolves to something.
    stages["main"] += hpccm.primitives.shell(
        commands=[
            "test -x /usr/bin/python || "
            "update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && "
            "/usr/bin/python --version"
        ]
    )

    # Note that the list of stages should be sorted in dependency order.
    for build_stage in stages.values():
        if build_stage is not None:
            yield build_stage


if __name__ == "__main__":
    args = parser.parse_args()
    import hpccm.config

    # Set container specification output format
    hpccm.config.set_container_format(args.format)
    # Normally in hpccm the first call to baseimage sets the context
    # for other packages, e.g. which distro version to use. We want to
    # set that early on, so that hpccm can do the right thing to use
    # the right llvm apt repo when we want to use versions of llvm
    # that were never supported by the main apt repo.
    hpccm.config.set_linux_distro(hpccm_distro_name(args))

    container_recipe = build_stages(args)

    # Output container specification
    for stage in container_recipe:
        print(stage)
