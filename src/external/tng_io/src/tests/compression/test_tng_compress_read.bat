@echo off
setlocal enableextensions enabledelayedexpansion
SET /A I=0
:start
SET /A I+=1
test_tng_compress_read%I%
IF "%I%" == "78" (
  GOTO end
) ELSE (
  GOTO start
)
:end
endlocal
