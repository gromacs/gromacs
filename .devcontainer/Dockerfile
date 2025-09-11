# ARG VARIANT=ubuntu-22.04
# FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

# FROM registry.gitlab.com/gromacs/gromacs/ci-ubuntu-20.04-llvm-9-cuda-11.0.3
# FROM registry.gitlab.com/gromacs/gromacs/ci-ubuntu-20.04-llvm-9-docs

FROM registry.gitlab.com/gromacs/gromacs/ci-ubuntu-24.04-llvm-18-cuda-12.6.3-adaptivecpp-24.10.0-rocm-6.3.1

LABEL maintainer="Bernd Doser <bernd.doser@h-its.org>"

# ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="3.31.3"

# # Optionally install the cmake for vcpkg
# COPY ./reinstall-cmake.sh /tmp/

# RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
#         chmod +x /tmp/reinstall-cmake.sh && /tmp/reinstall-cmake.sh ${REINSTALL_CMAKE_VERSION_FROM_SOURCE}; \
#     fi \
#     && rm -f /tmp/reinstall-cmake.sh

RUN apt-get update \
 && export DEBIAN_FRONTEND=noninteractive \
 && apt-get -y install --no-install-recommends \
       doxygen \
      #  libomp-14-dev \
       python3-pip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# RUN pip install --upgrade pip \
#  && hash -r pip3

# RUN pip install \
#       dash \
#       pandas

ARG USERNAME=vscode
ARG USER_UID=15005
ARG USER_GID=15000

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# ********************************************************
# * Anything else you want to do like clean up goes here *
# ********************************************************

# [Optional] Set the default user. Omit if you want to keep the default as root.
USER $USERNAME