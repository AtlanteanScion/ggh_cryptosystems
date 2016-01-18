@echo off

md obj
cd obj

:: Look for includes in ../include/
:: Compile; don't link
g++ -I ../include -c ../src/*.cpp

cd ..

:: The archiver can be used to place object files into a static library.
:: r = Insert object files into archive (with replacement).
:: c = Create the archive, if it doesn't already exist.
:: s = Maintain the table that maps symbol names to object file names.
ar rcs libntl.a obj/*.o

pause

