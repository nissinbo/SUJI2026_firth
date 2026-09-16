param([ValidateSet('pilot','full')][string]$Phase='full',[ValidateSet('proc','fl')][string]$Engine='proc',[int]$From=1,[int]$To=420,[string]$SasPath="C:/Program Files/SASHome/SASFoundation/9.4/sas.exe")
$ErrorActionPreference='Stop'
$analysisRoot=(Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$analysisStem='sas_'+$Phase+'_'+$Engine
$analysisOutput=Join-Path $analysisRoot ('results/raw/'+$analysisStem)
if ($Engine -eq 'fl') { $analysisStem += '_'+$From+'_'+$To }
$analysisWork=Join-Path $analysisRoot ('tmp/runtime/'+$analysisStem+'_work')
New-Item -ItemType Directory -Force -Path $analysisOutput,$analysisWork | Out-Null
$env:FIRTH_PROJECT_ROOT=$analysisRoot.Replace('\','/')
$env:FIRTH_PHASE=$Phase
$env:FIRTH_ENGINE=$Engine
$env:FIRTH_FROM=[string]$From
$env:FIRTH_TO=[string]$To
$analysisSource=if ($Engine -eq 'fl') { 'sas/fit_fl.sas' } else { 'sas/fit_logistic.sas' }
$analysisArgs=@('-sysin',('"'+(Join-Path $analysisRoot $analysisSource)+'"'),'-log',('"'+(Join-Path $analysisRoot ('tmp/runtime/'+$analysisStem+'.log'))+'"'),'-print',('"'+(Join-Path $analysisRoot ('tmp/runtime/'+$analysisStem+'.lst'))+'"'),'-work',('"'+$analysisWork+'"'),'-nosplash','-noterminal','-nologo')
$analysisProcess=Start-Process -FilePath $SasPath -ArgumentList $analysisArgs -WindowStyle Hidden -PassThru
$analysisProcess.Id | Set-Content -LiteralPath (Join-Path $analysisRoot ('tmp/runtime/'+$analysisStem+'.pid'))
Write-Output "Started $analysisStem PID=$($analysisProcess.Id)"
$analysisProcess.WaitForExit()
Write-Output "Completed $analysisStem exit=$($analysisProcess.ExitCode)"
if ($analysisProcess.ExitCode -gt 1) { throw "SAS failed with exit code $($analysisProcess.ExitCode)" }
