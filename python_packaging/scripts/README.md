# Testing scripts

The scripts in this directory are for convenience when running different suites of tests.
They can be specified as arguments to the Docker containers based on `ci.dockerfile`

* `integrationtest` runs the test suite in the higher level `test` directory
* `required_functionality` runs the script extracted from `acceptance/RequiredFunctionality.ipynb`
* `run_pytest` runs single-threaded tests for the gmxapi Python package.
* `run_pytest_mpi` launches a 2-rank MPI session for the gmxapi Python package tests.
