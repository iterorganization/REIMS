# If EcosimPro is still open in the server, close it
# Get-Process
$processes = Get-Process
$process = $processes | Where-Object {$_.Name -eq "EcosimPro"} 
if($process){
	Write-Output "Existing EcosimPro Process Detected. Stopping Process."
    Stop-Process -Name "EcosimPro"
}  

# Mount the K drive
If (!(Test-Path K:)) {
	Write-Output "Mounting"
    net use K: \\iter.org\cfs\group\Plant_modelling\CI_store /Y /persistent:yes
}

# Ensure the K drive is mounted to the right path.
$currentPath = Get-PSDrive K
Write-Output $currentPath.DisplayRoot
If ($currentPath.DisplayRoot -ne "\\iter.org\cfs\group\Plant_modelling\CI_store") {
    net use K: \\iter.org\cfs\group\Plant_modelling\CI_store /Y /persistent:yes
}

Write-Output "K: Drive mounted"