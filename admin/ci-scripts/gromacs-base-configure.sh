#!/usr/bin/env bash
set -e
CMAKE=${CMAKE:-$(which cmake)}
echo "Running on host:" $KUBERNETES_HOSTNAME
echo $CMAKE_COMPILER_SCRIPT
echo $CMAKE_EXTRA_OPTIONS
echo $CMAKE_SIMD_OPTIONS
echo $CMAKE_GPU_OPTIONS
echo $CMAKE_MPI_OPTIONS
echo $CMAKE_PRECISION_OPTIONS
echo $CMAKE_BUILD_TYPE_OPTIONS
echo $CMAKE_GMXAPI_OPTIONS
if [[ -d $BUILD_DIR ]] ; then
      rm -rf $BUILD_DIR && mkdir $BUILD_DIR ;
else
      echo "Preparing new build directory" ;
      mkdir $BUILD_DIR
fi
cd $BUILD_DIR
which $CMAKE
$CMAKE --version
set +e -o pipefail # Make "$?" work correctly
$CMAKE .. \
      -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
      $CMAKE_COMPILER_SCRIPT \
      $CMAKE_EXTRA_OPTIONS \
      $CMAKE_SIMD_OPTIONS \
      $CMAKE_MPI_OPTIONS \
      $CMAKE_PRECISION_OPTIONS \
      $CMAKE_BUILD_TYPE_OPTIONS \
      $CMAKE_GPU_OPTIONS \
      $CMAKE_GMXAPI_OPTIONS \
      -DCMAKE_INSTALL_PREFIX=../$INSTALL_DIR -DGMX_COMPILER_WARNINGS=ON \
      2>&1 | tee cmakeLog.log

EXITCODE=$?
set -e +o pipefail

if [[ "$CMAKE_GPU_OPTIONS" == *"GMX_SYCL=ACPP"* ]]; then
    # The GROMACS CMake code sets some non-cache variables
    # that are then set as non-forced cache variables by the
    # AdaptiveCPP CMake config file, which seems to be
    # recognized as two variables of the same name by
    # _getACppCmakeFlags() during a non-initial run of CMake.
    # Those copies are then added to _ALL_ACPP_CMAKE_FLAGS
    # which are then cached with different values from the
    # first run, some of which are duplicates.
    #
    # Maybe this can be improved but should not be attempted
    # until after support for hipSYCL is dropped. For now, the
    # expectation that a second run of CMake has a stable
    # CMakeCache.txt is not applied to AdaptiveCPP builds.
    #
    # See issues #4716 and #4720 for related content.
    :
elif [[ "$CMAKE_EXTRA_OPTIONS" == *"GMX_NNPOT=TORCH"* ]]; then
    # The find package from Pytorch changes which CUDA libraries it
    # uses between /usr/local/cuda and /usr/local/cuda/x.y between the
    # first and second call of CMake, so we can't require that the
    # CMake cache is stable in this case. Its CMake code is also
    # unilaterally talkative, so also we don't require that a second
    # call of CMake is quiet in this case.
    :
elif [[ "$CMAKE_EXTRA_OPTIONS" == *"GMX_GPU_FFT_LIBRARY=VkFFT"* ]] && [[ "$CMAKE_EXTRA_OPTIONS" == *"HIP_TARGET_ARCH"* ]]; then
    # The use of VkFFT with a HIP target requires finding hiprtc
    # and the HIP-finding machinery we use is not robust enough
    # to handle this use case with a stable CMakeLists.txt while
    # handling VkFFT+ACPP cases correctly. So we don't require a
    # stable CMakeLists.txt in this case.
    :
elif [[ $EXITCODE ]] ; then
    # CMake completed successfully the first time. Prepare to check
    # that CMakeCache.txt does not change when CMake runs again.
    cp CMakeCache.txt initial-CMakeCache.txt
    touch CMakeCache.txt

    echo "Running CMake a second time after touching CMakeCache.txt..." >> cmakeLog.log
    $CMAKE .. 2>&1 | tee secondCMakeLog.log

    EXITCODE=$?
    if [[ $EXITCODE ]] ; then
        # CMake completed successfully the second time. Now we check
        # the behavior. Note that e.g. the terminal output can be
        # affected by the use of e.g. CMAKE_WARN_DEPRECATED, which may
        # have been set for particular build configurations.

        echo "Checking that CMakeCache.txt did not change when run a second time. Diff was:" >> cmakeLog.log
        # Quiet diff just to see whether there is a difference. If there
        # is one, we report it later.
        diff initial-CMakeCache.txt CMakeCache.txt
        EXITCODE=$?
        # When there was a difference, diff exits with a code of 1
        if [[ "$EXITCODE" -eq 0 ]] ; then
            echo "empty! That's good!"
            EXITCODE=0
        else
            # Take the diff again to capture the diff to the log file
            diff initial-CMakeCache.txt CMakeCache.txt | tee -a cmakeLog.log
            echo "The second run of cmake changed CMakeCache.txt, and it should not. Check the diff and fix the problem" | tee -a cmakeErrors.log
        fi

        # When suitable, check that the subsequent call to CMake was
        # suitably quiet on the terminal.
        if [[ "$CMAKE_GPU_OPTIONS" == *"GMX_GPU=HIP"* ]]; then
            # The ROCm CMake config files report status messages
            # unilaterally, so we don't require that a second call of
            # CMake is quiet in that case.
            :
        else
            # Do the check
            line_count=$(wc -l < "secondCMakeLog.log")
            if [ "$line_count" -ne 3 ] ; then
                echo "The log from the second run of CMake had extra output after the first three lines." >> cmakeErrors.log
                echo "It's OK to write status output during the first run of CMake, but thereafter it should be quiet" >> cmakeErrors.log
                echo "It contained:" >> cmakeErrors.log
                cat secondCMakeLog.log | tee -a cmakeErrors.log
            fi
        fi
    fi
fi

awk '/CMake Warning/,/^--|^$/' cmakeLog.log | tee -a cmakeErrors.log
awk '/CMake Error/,/^--|^$/' cmakeLog.log | tee -a cmakeErrors.log
if [ -s cmakeErrors.log  ] || [ $EXITCODE != 0 ]; then echo "Found CMake warning or error while processing build"; cat cmakeErrors.log ; exit 1; fi
cd ..
