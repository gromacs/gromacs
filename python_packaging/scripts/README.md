# Testing scripts

The scripts in this directory are for convenience when running different suites of tests.
They can be specified as arguments to the Docker containers based on `ci.dockerfile`

* `run_flake8` runs the Flake8 Python linting tool on the `gmxapi` package sources.
* `run_acceptace_test` runs the test suite in the higher level `test` directory.
* `run_gmxapi_unittest` runs single-threaded tests for the gmxapi Python package.
* `run_sample_test` runs the tests for the sample_restraint MD plugin API client code.
* `run_full` runs all of the above tests.
* `run_full_mpi` launches a 2-rank MPI session for the various Python tests.
