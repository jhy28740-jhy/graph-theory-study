$ErrorActionPreference = 'Continue'
$client = '178c6fc778ccc68e1d6a'
$scope  = 'repo%20read%3Aorg%20gist'
$tokFile = Join-Path $env:TEMP 'ghflow_token.txt'
if (Test-Path -LiteralPath $tokFile) { Remove-Item -LiteralPath $tokFile -Force }

$devJson = $null
for ($attempt = 1; $attempt -le 3; $attempt++) {
  $r = curl.exe -s --max-time 30 --retry 15 --retry-all-errors --retry-delay 2 --retry-max-time 300 -X POST 'https://github.com/login/device/code' -d ("client_id=" + $client + "&scope=" + $scope) -H 'Accept: application/json' 2>$null
  if ($r -match '"device_code"') { $devJson = $r | ConvertFrom-Json; break }
  Start-Sleep -Seconds 3
}
if (-not $devJson) { Write-Output 'FLOW-FAIL: cannot obtain device code'; exit 1 }
Write-Output ("USER-CODE=" + $devJson.user_code)
Write-Output ("VERIFY=" + $devJson.verification_uri)

$deadline = (Get-Date).AddMinutes(14)
$token = $null
while ((Get-Date) -lt $deadline) {
  $r = curl.exe -s --max-time 30 --retry 12 --retry-all-errors --retry-delay 2 --retry-max-time 240 -X POST 'https://github.com/login/oauth/access_token' -d ("client_id=" + $client + "&device_code=" + $devJson.device_code + "&grant_type=urn:ietf:params:oauth:grant-type:device_code") -H 'Accept: application/json' 2>$null
  if ($r -match '"access_token"') { $token = ($r | ConvertFrom-Json).access_token; break }
  if ($r -match '"(expired_token|access_denied|unsupported_grant_type|incorrect_device_code)"') { Write-Output ('FLOW-FAIL: ' + $r); exit 1 }
  if ($r -match '"slow_down"') { Start-Sleep -Seconds 10 } else { Start-Sleep -Seconds 5 }
}
if (-not $token) { Write-Output 'FLOW-FAIL: timeout waiting for authorization'; exit 1 }
Set-Content -LiteralPath $tokFile -Value $token -NoNewline
Write-Output 'FLOW-OK'
exit 0
