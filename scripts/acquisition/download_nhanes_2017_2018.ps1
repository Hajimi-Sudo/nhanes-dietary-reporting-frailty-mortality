$ErrorActionPreference = 'Stop'

$projectRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$downloadRoot = Join-Path $projectRoot 'data\raw\nhanes_2017_2018'
$logRoot = Join-Path $projectRoot 'logs'
New-Item -ItemType Directory -Path $downloadRoot, $logRoot -Force | Out-Null

$downloadDate = Get-Date -Format 'yyyy-MM-dd'
$headers = @{ UserAgent = 'survey-frailty-reproducibility/0.1' }
$nhanesBase = 'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles'
$mortalityBase = 'https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality'
$files = @(
    @{ Name = 'DEMO_J.xpt'; Url = "$nhanesBase/DEMO_J.xpt" },
    @{ Name = 'DR1IFF_J.xpt'; Url = "$nhanesBase/DR1IFF_J.xpt" },
    @{ Name = 'DR1TOT_J.xpt'; Url = "$nhanesBase/DR1TOT_J.xpt" },
    @{ Name = 'BMX_J.xpt'; Url = "$nhanesBase/BMX_J.xpt" },
    @{ Name = 'PFQ_J.xpt'; Url = "$nhanesBase/PFQ_J.xpt" },
    @{ Name = 'MCQ_J.xpt'; Url = "$nhanesBase/MCQ_J.xpt" },
    @{ Name = 'DPQ_J.xpt'; Url = "$nhanesBase/DPQ_J.xpt" },
    @{ Name = 'BPQ_J.xpt'; Url = "$nhanesBase/BPQ_J.xpt" },
    @{ Name = 'DIQ_J.xpt'; Url = "$nhanesBase/DIQ_J.xpt" },
    @{ Name = 'HUQ_J.xpt'; Url = "$nhanesBase/HUQ_J.xpt" },
    @{ Name = 'RXQ_RX_J.xpt'; Url = "$nhanesBase/RXQ_RX_J.xpt" },
    @{ Name = 'KIQ_U_J.xpt'; Url = "$nhanesBase/KIQ_U_J.xpt" },
    @{ Name = 'GHB_J.xpt'; Url = "$nhanesBase/GHB_J.xpt" },
    @{ Name = 'CBC_J.xpt'; Url = "$nhanesBase/CBC_J.xpt" },
    @{ Name = 'PAQ_J.xpt'; Url = "$nhanesBase/PAQ_J.xpt" },
    @{ Name = 'SMQ_J.xpt'; Url = "$nhanesBase/SMQ_J.xpt" },
    @{ Name = 'ALQ_J.xpt'; Url = "$nhanesBase/ALQ_J.xpt" },
    @{ Name = 'NHANES_2017_2018_MORT_2019_PUBLIC.dat'; Url = "$mortalityBase/NHANES_2017_2018_MORT_2019_PUBLIC.dat" }
)

$manifest = foreach ($file in $files) {
    $target = Join-Path $downloadRoot $file.Name
    Invoke-WebRequest -Uri $file.Url -Headers $headers -UseBasicParsing -OutFile $target
    $item = Get-Item -LiteralPath $target
    if ($item.Length -lt 1000) { throw "Downloaded file is unexpectedly small: $target" }
    $prefix = [System.Text.Encoding]::ASCII.GetString([System.IO.File]::ReadAllBytes($target)[0..7])
    if ($prefix -match '<!DOCTYPE|<html') { throw "Downloaded HTML instead of data: $target" }
    $hash = Get-FileHash -LiteralPath $target -Algorithm SHA256
    [pscustomobject]@{
        dataset = 'NHANES 2017-2018 and linked mortality'
        file = $file.Name
        source_url = $file.Url
        downloaded_at = $downloadDate
        bytes = $item.Length
        sha256 = $hash.Hash
    }
}

$manifestPath = Join-Path $logRoot 'nhanes_download_manifest.csv'
$manifest | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding UTF8
Write-Output "Downloaded $($manifest.Count) files"
Write-Output $manifestPath





