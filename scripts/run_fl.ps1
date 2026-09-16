$ErrorActionPreference='Stop'
New-Item -ItemType Directory -Force -Path (Join-Path $PSScriptRoot '../tmp/runtime') | Out-Null
$flParallelRoot=(Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$flParallelJobs=@()
foreach ($flRange in @(@(1,210),@(211,280),@(281,350),@(351,420))) {
 $flArgs=@('-ExecutionPolicy','Bypass','-File',('"'+(Join-Path $flParallelRoot 'scripts/run_sas.ps1')+'"'),'-Phase','full','-Engine','fl','-From',[string]$flRange[0],'-To',[string]$flRange[1])
 $flParallelJobs += Start-Process -FilePath 'powershell.exe' -ArgumentList $flArgs -WindowStyle Hidden -PassThru -WorkingDirectory $flParallelRoot -RedirectStandardOutput (Join-Path $flParallelRoot ('tmp/runtime/fl_runner_'+$flRange[0]+'.txt')) -RedirectStandardError (Join-Path $flParallelRoot ('tmp/runtime/fl_runner_'+$flRange[0]+'.err'))
}
Write-Output 'Started four disjoint FL ranges; saved results will be reused.'
foreach ($flJob in $flParallelJobs) { $flJob.WaitForExit() }
if (@($flParallelJobs | Where-Object {$_.ExitCode -ne 0}).Count -gt 0) { throw 'An FL range failed; inspect tmp/runtime/fl_runner_*.err.' }
Write-Output 'All FL ranges finished. Individual numerical failures remain recorded in results.'
