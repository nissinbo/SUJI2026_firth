$ErrorActionPreference='Stop'
New-Item -ItemType Directory -Force -Path (Join-Path $PSScriptRoot '../tmp/runtime') | Out-Null
$analysisRoot=(Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$analysisJobs=@()
foreach ($analysisShard in 1..2) {
 foreach ($analysisRange in @(@(1,53),@(54,106),@(107,158),@(159,210))) {
  $analysisLabel='r_'+$analysisShard+'_'+$analysisRange[0]+'_'+$analysisRange[1]
  $analysisRArgs=@('"'+(Join-Path $analysisRoot 'R/run_models.R')+'"', '--shard='+$analysisShard,'--from='+$analysisRange[0],'--to='+$analysisRange[1])
  $analysisJobs += Start-Process -FilePath (Get-Command Rscript -ErrorAction Stop).Source -ArgumentList $analysisRArgs -WindowStyle Hidden -PassThru -WorkingDirectory $analysisRoot -RedirectStandardOutput (Join-Path $analysisRoot ('tmp/runtime/'+$analysisLabel+'.log')) -RedirectStandardError (Join-Path $analysisRoot ('tmp/runtime/'+$analysisLabel+'.err'))
 }
}
$analysisJobs.Id | Set-Content -LiteralPath (Join-Path $analysisRoot 'tmp/runtime/r_worker_pids.txt')
Write-Output "Started $($analysisJobs.Count) R workers using disjoint checkpoint ranges."
foreach ($analysisJob in $analysisJobs) { $analysisJob.WaitForExit() }
if (@($analysisJobs | Where-Object {$_.ExitCode -ne 0}).Count -gt 0) { throw 'An R worker failed; inspect tmp/runtime.' }
Write-Output 'All R workers completed.'
