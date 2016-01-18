@echo off

:: Look for includes in include/
:: Look for libraries in current directory
:: Use library libntl.a
g++ -I include -L . %1 -l ntl -std=c++11 -o %~n1.exe

pause

