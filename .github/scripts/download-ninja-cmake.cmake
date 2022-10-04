set(cmake_version $ENV{CMAKE_VERSION})
set(ninja_version $ENV{NINJA_VERSION})

message(STATUS "Using host CMake version: ${CMAKE_VERSION}")
message(STATUS "Using RUNNER_OS: $ENV{RUNNER_OS}")

if ("$ENV{RUNNER_OS}" STREQUAL "Windows")
    set(ninja_suffix "win.zip")
    set(cmake_suffix "win64-x64.zip")
    set(cmake_dir "cmake-${cmake_version}-win64-x64/bin")
elseif ("$ENV{RUNNER_OS}" STREQUAL "macOS")
    set(ninja_suffix "mac.zip")
    set(cmake_suffix "Darwin-x86_64.tar.gz")
    set(cmake_dir "cmake-${cmake_version}-Darwin-x86_64/CMake.app/Contents/bin")
elseif ("$ENV{RUNNER_OS}" STREQUAL "Linux" AND "$ENV{RUNNER_ARCH}" STREQUAL "ARM64")
    set(ninja_suffix "linux.zip")
    set(cmake_suffix "linux-aarch64.tar.gz")
    set(cmake_dir "cmake-${cmake_version}-linux-aarch64/bin")
endif()

set(ninja_url "https://github.com/ninja-build/ninja/releases/download/v${ninja_version}/ninja-${ninja_suffix}")
file(DOWNLOAD "${ninja_url}" ./ninja.zip SHOW_PROGRESS)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ./ninja.zip)

set(cmake_url "https://github.com/Kitware/CMake/releases/download/v${cmake_version}/cmake-${cmake_version}-${cmake_suffix}")
file(DOWNLOAD "${cmake_url}" ./cmake.zip SHOW_PROGRESS)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ./cmake.zip)

# Add to PATH environment variable
file(TO_CMAKE_PATH "$ENV{GITHUB_WORKSPACE}/${cmake_dir}" cmake_dir)
set(path_separator ":")
if ("$ENV{RUNNER_OS}" STREQUAL "Windows")
    set(path_separator ";")
endif()
file(APPEND "$ENV{GITHUB_PATH}" "$ENV{GITHUB_WORKSPACE}${path_separator}${cmake_dir}")

if (NOT "$ENV{RUNNER_OS}" STREQUAL "Windows")
    execute_process(
    COMMAND chmod +x ninja
    COMMAND chmod +x ${cmake_dir}/cmake
    )
endif()
