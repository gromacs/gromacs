# Acceptance testing

`../roadmap.rst` summarizes the expected functionality to deliver during 2019 Q2.
`RequiredFunctionality.ipynb` focuses on functionality and behavior related 
to the Python package.

Targeted behavior can be explored and verified interactively as a
[Jupyter notebook](https://jupyter.org)
or a runnable Python script can be extracted.

## Accessing the notebook

### Download a Docker image

For a feature under review, we can push images to DockerHub for the convenience
of reviewers.

For a feature request tagged as `fr1`, the following `docker` command line
pulls the `acceptance` repository tagged `fr1` from the `gmxapi` project
namespace on
[dockerhub](Run a container and map its notebook server to a local port.)
(if needed and available)
and maps port 8888 of the container to port 8888 on `localhost`.
Without additional arguments, the container will run a jupypter notebook server
and print to the terminal a URL that can be browsed to on the machine where
`docker` is run.

    docker run --rm -t -p 8888:8888 gmxapi/acceptance:fr1

Browse to `RequiredFunctionality.ipynb` in the `acceptance` directory

Note that changes made in this session will not be retained unless you use the
web interface to export a new `ipynb` to your machine
(such as to replace the repository version).
Please clear the output and/or use the `strip_notebook.py` script to trim
the notebook before committing back to the git repository.

### Build a Docker image

Use `../docker/acceptance.dockerfile` image as documented in
`../docker/README.md` and in the Dockerfile itself.

### Build it all yourself

1. Build and install GROMACS with GMXAPI=ON
2. Source the GMXRC
3. Create and activate a Python virtualenv. (optional but recommended)
4. `pip install` the `gmxapi` Python package from `python_packaging/gmxapi`
5. Run a Jupyter notebook server and navigate to RequiredFunctionality.ipynb in this directory.

### GitHub

GitHub nicely displays the notebook at (for example)
https://github.com/eirrgang/gromacs-gmxapi/blob/fr2/python_packaging/acceptance/RequiredFunctionality.ipynb

### Convert to some other format

Install Jupyter (e.g. `pip install jupyter`) and use its `nbconvert` utility.

Examples:

    jupyter nbconvert RequiredFunctionality.ipynb --to python
    jupyter nbconvert RequiredFunctionality.ipynb --to pdf

See `jupyter nbconvert --help` for details.
