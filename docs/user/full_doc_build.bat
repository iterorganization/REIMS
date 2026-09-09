@echo off
if not defined VIRTUAL_ENV (
    call activate global_python
)
rem Add MiKTexX and perl to the PATH if not exist yet
set PATH=%PATH%;C:\Soft\MiKTeX\miktex\bin\x64;C:\Soft\perl-5.38.2.2\perl\bin

rem This script builds the full documentation for the project.
cd ..\..\tools
python yaml_to_json.py
cd ..\docs\user
echo Building PDF documentation (user manual)...
call make.bat latex
python fix_latex.py
cd _build\latex
pdflatex reims.tex
pdflatex reims.tex
cd ..\..

echo Building HTML documentation (user manual)...
call make.bat html

echo Building DLL developer guide (HTML)...
cd ..\dll
call make.bat html

echo Building DLL developer guide (PDF)...
call make.bat latex
cd _build\latex
pdflatex reimsdlldeveloperguide.tex
pdflatex reimsdlldeveloperguide.tex
cd ..\..\..
cd user
