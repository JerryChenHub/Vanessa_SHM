param([int]$SolidWorksPid=8704)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path $PSScriptRoot -Parent
$cadDir=Join-Path $projectRoot 'models\skin_demo\cad'
New-Item -ItemType Directory -Path $cadDir -Force | Out-Null
$hashes=@()
foreach($tag in @('1','2')) {
    $source=Join-Path 'C:\Users\Jerry12138\Downloads' ($tag+'.SLDPRT')
    $dest=Join-Path $cadDir ($tag+'.SLDPRT')
    Copy-Item -LiteralPath $source -Destination $dest
    $hashes+=@{source=$source;source_sha256=(Get-FileHash -LiteralPath $source).Hash;copy=$dest;copy_sha256=(Get-FileHash -LiteralPath $dest).Hash}
}
$interop='D:\sw2019\SOLIDWORKS\api\redist\SolidWorks.Interop.sldworks.dll'
Add-Type -Path $interop
Add-Type -ReferencedAssemblies @($interop,'System.Core') -Path (Join-Path $PSScriptRoot 'ExportSkinGeometry.cs')
[ExportSkinGeometry]::Run($SolidWorksPid,$cadDir) | ConvertTo-Json -Depth 10 | Set-Content -Encoding UTF8 (Join-Path $cadDir 'export.json')
foreach($h in $hashes) {
    if((Get-FileHash -LiteralPath $h.source).Hash -ne $h.source_sha256){throw 'Source CAD changed'}
}
$hashes | ConvertTo-Json -Depth 5 | Set-Content -Encoding UTF8 (Join-Path $cadDir 'provenance.json')
Get-ChildItem -LiteralPath $cadDir | Select-Object Name,Length
