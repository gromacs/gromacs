@echo off
set DIRNAME=%~dp0%
REM BAT files only see args 1-9. We use shift to get all.
set ARGS=
:loop
if [%1] == [] goto endloop
set ARGS=%ARGS% %1
shift
goto loop
:endloop
%DIRNAME%gmx-ana.exe -type sas %ARGS%
