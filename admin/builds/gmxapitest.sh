#!/bin/bash

#Read in version as run time variable
version=$1
#Read in C compiler as runtime variable
cc=$2
#Read in CXX compiler as runtime variable
cxx=$3
#Read in directory that contains the GMXRC script
gmxrc="$4"
echo "gmxapi version is ${version}"
echo "C compiler is $cc"
echo "CXX compiler is $cxx"
echo "Location of GMXRC script is $gmxrc"

export CC=$cc
export CXX=$cxx

# Checks return code and cats output files to terminal if
# the tests failed
checkReturnCode() {

returnCode=$1
file=$2

if [[ "$returnCode" == "1" ]]; then
    echo "Return code was $returnCode"
    echo "Content of $file"
    cat ./$file
fi
}

# Run the actual test of gmxapi in the chosen virtual environment
# Argument is the python venv, script returns early with return code 4
# if no argument is given (should not happen in Jenkins). Other return codes
# are selected according to the result of the tests themselves.
function runTest() {

pythonVenv=$1
echo "Python environment is $pythonVenv"
if [ -z $pythonVenv ] ; then
    echo "No argument given to launch venv"
    return 4
fi

if [[ "$pythonVenv" != "py2env" ]] && [[ "$pythonVenv" != "py3env" ]] ; then
    echo "Unknown virtualenvironment"
    return 5
fi

venvLine="$HOME/${pythonVenv}/bin/activate"
source $venvLine
export PYTHON=`which python`
echo "Running under $PYTHON"

# Prepare gmxapi and plugin under venv with correct version handle
pip uninstall -y gmx #-$version
pip uninstall -y myplugin #-$version

# get contents of gmxapi directory to current wd
cp -r $HOME/gmxapi/gmxapi $HOME/gmxapi/sample_restraint .

wd=`pwd`
echo "Working in $wd"

(cd gmxapi && ./ci_scripts/pygmx.sh) &> pygmx.out
pygmx_test=$?
checkReturnCode $pygmx_test "pygmx.out"

(cd sample_restraint && ./ci_scripts/sample_restraint.sh) &> sample_restraint.out
sample_restraint_test=$?
checkReturnCode $sample_restraint_test "sample_restraint.out"
deactivate
return_code=$((pygmx_test+2*sample_restraint_test))

return $return_code

}

source $gmxrc

runTest "py2env"

returnCode=$?

if [[ "$returnCode" != "0" ]]; then
    exit $returnCode
fi

runTest "py3env"

returnCode=$?

if [[ "$returnCode" != "0" ]] ; then
    exit $returnCode
fi

