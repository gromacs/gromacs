## pybind11

For simplicity, the pybind headers were retrieved from 
https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz

    git rm -rf pybind
    mkdir -p pybind
    wget https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz
    tar xvf *.tar.gz
    mv pybind11*/include pybind/
    mv pybind11*/tools pybind/
    cp pybind11*/{CMakeLists.txt,LICENSE} pybind/
    git add pybind

If we need to update the version of pybind, the easiest thing would be to just `git rm` the pybind directory, 
decompress a new one, and `git add` it back (in a single commit).

Alternatives include git submodules and `git subtree`
(see https://www.atlassian.com/blog/git/alternatives-to-git-submodule-git-subtree), or to rely on local pybind 
sournces, possibly retrieved as a `pip` dependency.

### source tree pruning

To minimize the footprint of pybind, we include a minimal amount of it.

`docs` and `tests` directories removed
