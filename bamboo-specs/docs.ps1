#$out_path = 'ci_test'
#$doc_path = 'ford_generated'

$lastExitCode = 0

# Setting shorter variable names from the bamboo env requirements
$7z = $env:bamboo_capability_system_builder_command_7Zip

Write-host "Mounting K: from '\\iter.org\cfs\group\Plant_modelling\CI_store'"
If (!(Test-Path k:)) {
    net use k: \\iter.org\cfs\group\Plant_modelling\CI_store
    If ($lastExitCode -ne 0){throw "Can't mount network drive"}
}

# Setting new path 
$path  = "k:\Tools\graphviz;"
$path += "k:\Tools\ford;"
$path += "k:\Tools\ford\Library\mingw-w64\bin;"
$path += "k:\Tools\ford\Library\usr\bin;"
$path += "k:\Tools\ford\Library\bin;"
$path += "k:\Tools\ford\Scripts;"
$env:path = "$path$env:path"

Write-host 'Running FORD'
Set-Location docs
python k:/Tools/ford/Scripts/ford.exe ford_doc.md
If ($lastExitCode -ne 0){throw "Ford exit with errors"}

# Artefact generation
& $7z a -r -y ford.zip ford_generated\*.* > out.txt

# ----------------------------------------------------------------
#             This should be on deployment server
# ----------------------------------------------------------------

# url: https://sharepoint.iter.org/units/as/system%20engineering
$NetPath  = "\\sharepoint.iter.org@ssl\units\as\system engineering"
net use $NetPath

If (Test-Path "$NetPath\reims\ford_generated"){
    Try{
        Remove-item "$NetPath\reims\ford_generated" -Force -Recurse
        Write-Host "Previous files are removed"
    }
    Catch [Exception]{}
}

Write-host "Start copy of new html files"
copy-item  ford_generated -Destination $NetPath\reims -Recurse
Write-host "Copy done"
