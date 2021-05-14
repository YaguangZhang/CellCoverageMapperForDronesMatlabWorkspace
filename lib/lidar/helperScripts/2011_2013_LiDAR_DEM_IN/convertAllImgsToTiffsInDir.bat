:: Convert all the .img file to .tif (recursively) in the input dir.
::     - The input dir is default to ".";
::     - This script has to be run in the OSGeo4W shell.

@echo off
SET startpath=%1

IF "%startpath%"=="" (
  SET startpath="."
  )

ECHO Target director is %startpath%
PAUSE

FOR /R %startpath% %%f IN (*.img) do (call :translate "%%f")
PAUSE
GOTO :eof

:translate
  SET var=%1
  SET curfile=%var% 
  SET targetfile=%var:~0,-5%.tif"
  ECHO %curfile%
  IF EXIST %targetfile%.aux.xml (
    ECHO Already translated.
    ECHO Skipped.
  ) ELSE (
    ECHO Translating ...
    gdal_translate %curfile% %targetfile%
    ECHO Done!
  )  