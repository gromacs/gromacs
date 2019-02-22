# Acceptance testing

`../roadmap.rst` summarizes the expected functionality to deliver during 2019 Q2. `RequiredFunctionality
.ipynb` focuses on functionality and behavior related 
to the Python package.

Targeted behavior can be explored and verified interactively.

## Accessing the notebook

### Download a Docker image

For a feature under review, there will be images available on DockerHub.
Run a container and map its notebook server to a local port.
Example:

    docker run --rm -t -p 8888:8888 gmxapi/acceptance:fr1

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
