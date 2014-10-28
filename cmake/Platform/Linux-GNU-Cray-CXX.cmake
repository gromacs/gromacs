include(Platform/Linux-GNU-Cray)
__linux_compiler_gnu(C)
__linux_compiler_gnu(CXX)

# Set up cache variables that should always be this way (no need to
# use FORCE, the user might know what they are doing).
set(CMAKE_SKIP_RPATH ON CACHE BOOL "Can't use RPATH with Cray")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Building statically on Cray")
set(GMX_PREFER_STATIC_LIBS ON CACHE BOOL "Building statically on Cray")
set(GMX_LINK_STATIC_BINARIES ON CACHE BOOL "Avoid linking shared versions of e.g. system libraries")

