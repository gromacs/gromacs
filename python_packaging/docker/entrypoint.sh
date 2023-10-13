#!/usr/bin/env bash
set -e

# Set up user environment here.
cd $HOME
export PATH=/docker_entry_points:$PATH
# Activate the Python venv, if present for the current user.
# (testing user is default, but terminal user can be specified to `docker run`)
# Note: VENV defined in ci.dockerfile
test -f $VENV/bin/activate && source $VENV/bin/activate

if [ "$1" -a ! -x "$1" -a ! "$(which $1)" ]; then
    echo "Cannot execute $1"
    echo "Use an installed executable or one of the entry point scripts:"
    ls -1 /docker_entry_points | grep 'run_'
    exit 1
fi
exec "$@"
