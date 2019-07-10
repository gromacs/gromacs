# Build a C++ library using libgmxapi to extend GROMACS.
# We process this subdirectory first to reaffirm that this code does not depend
# on other substantial infrastructure, such as a particular Python bindings
# framework or the Googletest framework.
add_subdirectory(cpp)

# Build a Python extension package from our new library.
add_subdirectory(pythonmodule)
