REM Requires Windows SDK and CMake to be in the path
REM Run Windows SDK SetEnv script
REM OpenMP is only available in SDK 6.1. Both SDKs support 32bit and 64bit.
REM Arbirtary choice to use 32bit with 6.1. 6.1 requires delayed expansion.
setlocal enabledelayedexpansion 
if "%GMX_OPENMP%"=="yes" (
   C:\Program Files\Microsoft SDKs\Windows\v6.1\Bin\SetEnv /Release /x86 
) else (
   C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv /Release
)
cmake -G "Visual Studio 10 Win64" -D GMX_DOUBLE=%GMX_DOUBLE% -D GMX_MPI=%GMX_MPI% -D GMX_OPENMP=%GMX_OPENMP% -DGMX_DEFAULT_SUFFIX=off . && ^
msbuild /m:2 /p:Configuration=MinSizeRel All_Build.vcxproj && ^
ctest -D ExperimentalTest -C MinSizeRel -V
