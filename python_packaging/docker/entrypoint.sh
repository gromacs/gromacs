#!/bin/bash
set -e

# Set up environment for "testing" user here.
cd $HOME
export PATH=$HOME/scripts:$PATH
# Activate the Python venv, if present for the current user.
# (testing user is default, but terminal user can be specified to `docker run`)
test -f $HOME/testing/bin/activate && source $HOME/testing/bin/activate

exec "$@"
