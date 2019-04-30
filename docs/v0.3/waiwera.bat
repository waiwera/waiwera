::
:: https://www.dostips.com/DtTipsStringManipulation.php#Snippets.MapLookup
:: https://stackoverflow.com/questions/3973824/windows-bat-file-optional-argument-parsing#comment4250214_3973824

:: print help if no argument
@ECHO OFF
IF "%1"=="" GOTO :help
IF "%1"=="help" GOTO :help
IF "%1"=="--help" GOTO :help
IF "%1"=="/help" GOTO :help
IF "%1"=="/h" GOTO :help
IF "%1"=="-h" GOTO :help

:: get current directory
SET v=%cd%

:: translate drive, small case, colon to slash, eg D: to /c/
SET drv=%v:~0,2%
SET map=A:-/a/;B:-/b/;C:-/c/;D:-/d/;E:-/e/;F:-/f/;G:-/g/;H:-/h/;I:-/i/;J:-/j/;K:-/k/;L:-/l/;M:-/m/;N:-/n/;O:-/o/;P:-/p/;Q:-/q/;R:-/r/;S:-/s/;T:-/t/;U:-/u/;V:-/v/;W:-/w/;X:-/x/;Y:-/y/;Z:-/z/
CALL SET drv=%%map:*%drv%-=%%
SET drv=%drv:;=&rem.%
:: ECHO %drv%

:: translate slashes to unix
SET _pth=%v:~3,1000%
SET pth=%_pth:\=/%
:: ECHO %pth%

SET np=%2
IF "%2"=="" (
    set np=1
)

:: put back together
set pp=%drv%%pth%
:: ECHO %pp% %np%
docker run -v %pp%:/data -w /data waiwera-phusion-debian mpiexec -np %np% /home/mpirun/waiwera/dist/waiwera %1
GOTO :end

:help

ECHO.
ECHO Run Waiwera at the current directory:
ECHO      waiwera input_file.json
ECHO   (or specify number of CPU, default is 1)
ECHO      waiwera input_file.json 4
ECHO.

:end

