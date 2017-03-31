This is a fork of the main Gromacs project in which interface, API, and extensibility issues are being investigated.
This README.md file supplants the main README file to avoid merge conflicts while providing convenient documentation to the BitBucket repository browser.

To use, you need to fetch the pybind11 git submodule, too.

    git clone git@bitbucket.org:kassonlab/gromacs_api.git
    cd gromacs_api
    git checkout trajectory
    git submodule init
    git submodule update src/external/pybind11