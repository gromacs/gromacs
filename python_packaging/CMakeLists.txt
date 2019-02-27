# Unimplemented.
#
# To be driven by a higher-level CMake configuration, a CMakeLists.txt file at this level would
# 1. configure: copy the contents of `src` to the build directory, potentially configuring version info
# 2. build: run (cd $CMAKE_CURRENT_BINARY_DIR/src && python setup.py sdist),
#    producing $CMAKE_CURRENT_BINARY_DIR/src/dist/gmxapi-0.1.0.dev0+fr7.tar.gz or similar
#    (add_custom_command(OUTPUT $sdist_name ...) where sdist_name is determined with
#    `python setup.py --full-name` and a choice of argument for `python setup.py sdist --formats`
# 3. install: copy the source distribution to $GROMACS_INSTALL_DIR/sdist/ or something
# 4. test:
#    0. copy src/test
#    1. require GROMACS library target to be present or importable
#    2. create and activate a Python venv in the build dir
#    3. `pip install <the_sdist_archive_file>` into the venv
#    4. run `pytest copy_of_src/test`
