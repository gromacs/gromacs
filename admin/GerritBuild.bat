REM Requires Windows SDK and CMake to be in the path
REM OpenMP is only available in SDK 6.1. Both SDKs support 32bit and 64bit.
REM Arbirtary choice to use 32bit with 6.1 (default)
REM 6.1 requires delayed expansion. If not enabled restarting cmd with /v
set A=B
if not !A!==B ( cmd /v /c admin\GerritBuild )
if "%GMX_OPENMP%"=="on" (
   set SetEnvCmd="C:\Program Files\Microsoft SDKs\Windows\v6.1\Bin\SetEnv" /Release
   set Target="NMake Makefiles"
) else (
   set SetEnvCmd="C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv" /Release
   set Target="Visual Studio 10 Win64"
)
%SetEnvCmd% && ^
cmake -G %Target% -D GMX_DOUBLE=%GMX_DOUBLE% -D GMX_MPI=%GMX_MPI% -D GMX_OPENMP=%GMX_OPENMP% -DGMX_DEFAULT_SUFFIX=off . && ^
msbuild /m:2 /p:Configuration=MinSizeRel All_Build.vcxproj && ^
ctest -D ExperimentalTest -C MinSizeRel -V
