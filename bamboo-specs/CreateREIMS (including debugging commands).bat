@ echo off
echo "-------- Creating reims.exe ---------"

:: Remove any old build folders
if exist "K:\reims\build" (
	rmdir K:\reims\build /s /q
)

:: Mount K drive if not already mounted
START /WAIT powershell -ExecutionPolicy Bypass -File "%~dp0MountDrive.ps1"

K:
cd K:\Tools\PortableIntel

powershell -Command "Unblock-File -Path K:\Tools\VisualStudioCompiler_14.43.34808_22621\Common7\IDE\CommonExtensions\Microsoft\TestWindow\Microsoft.VisualStudio.TestWindow.dll"

echo "------ Initializing oneAPI environment -------"
:: setvars-portable script sets all the necessary environment variables. Run it once per CMD session.
if "%SETVARS_COMPLETED%" NEQ "1" (
	@call setvars-portable
)

cd K:\reims

echo "Creating build folder."
mkdir ".\build"

cd build

::dir K:\Tools\PortableIntel\compiler\latest\bin

echo "-------- Running Cmake ---------"
set FC=K:/Tools/PortableIntel/compiler/latest/bin/ifort.exe 
set CC=K:/Tools/VisualStudioCompiler_14.43.34808_22621/VC/Tools/MSVC/14.43.34808/bin/Hostx64/x64/cl.exe
::set
cmake .. -G "NMake Makefiles" -D CMAKE_EXE_LINKER_FLAGS="/machine:x64" -D CMAKE_Fortran_COMPILER="K:/Tools/PortableIntel/compiler/latest/bin/ifort.exe" -D CMAKE_C_COMPILER="K:/Tools/VisualStudioCompiler_14.43.34808_22621/VC/Tools/MSVC/14.43.34808/bin/Hostx64/x64/cl.exe" 
:: -D CMAKE_GENERATOR="NMake Files"
:: -D CMAKE_ROOT="K:/Tools/VisualStudioCompiler_14.43.34808_22621/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/share/cmake-3.30" -D CMAKE_COMMAND="K:\Tools\VisualStudioCompiler_14.43.34808_22621\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin"
:: -D CMAKE_Fortran_COMPILER="K:/Tools/PortableIntel/compiler/latest/bin/ifort.exe" -D CMAKE_C_COMPILER="K:/Tools/VisualStudioCompiler_14.43.34808_22621/VC/Tools/MSVC/14.43.34808/bin/Hostx64/x64/cl.exe"
:: -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install

echo "----- Running Cmake Build ------"
::ninja 
cmake --build . --config Release

:: The fmifast.exe is the artifact Bamboo will recognize
copy "K:\reims\build\src\reims.exe" "K:\Artifacts\reims"
:: Copy additional dll's needed to run reims
copy "K:\Tools\VisualStudioCompiler_14.43.34808_22621\Windows Kits\10\bin\10.0.22621.0\x64\ucrt\ucrtbased.dll" "K:\Artifacts\reims"
copy "K:\Tools\VisualStudioCompiler_14.43.34808_22621\VC\Redist\MSVC\14.42.34433\onecore\debug_nonredist\x64\Microsoft.VC143.DebugCRT\vcruntime140d.dll" "K:\Artifacts\reims"
copy "K:\Tools\PortableIntel\compiler\latest\bin\libmmd.dll" "K:\Artifacts\reims"
copy "K:\Tools\PortableIntel\compiler\latest\bin\libiomp5md.dll" "K:\Artifacts\reims"
copy "K:\Tools\PortableIntel\compiler\latest\bin\libmmdd.dll" "K:\Artifacts\reims"
copy "K:\Tools\PortableIntel\compiler\latest\bin\libifcoremdd.dll" "K:\Artifacts\reims"
copy "K:\Tools\PortableIntel\compiler\latest\bin\libifcoremd.dll" "K:\Artifacts\reims"