$t = (Get-Content -LiteralPath (Join-Path $env:TEMP 'ghflow_token.txt') -Raw).Trim()
$bodyFile = Join-Path $env:TEMP 'ghflow_body.json'
Write-Output ("file exists: " + (Test-Path -LiteralPath $bodyFile) + " len=" + (Get-Item -LiteralPath $bodyFile).Length)
$hex = ([System.IO.File]::ReadAllBytes($bodyFile)[0..15] | ForEach-Object { $_.ToString('x2') }) -join ' '
Write-Output ("first bytes: " + $hex)
$out = curl.exe -s -S --max-time 30 --retry 6 --retry-all-errors -X POST 'https://api.github.com/user/repos' -H ("Authorization: Bearer " + $t) -H 'Accept: application/vnd.github+json' -H 'Content-Type: application/json' --data-binary ('@' + $bodyFile) -w "`nHTTP=%{http_code} UPLOAD=%{size_upload}" 2>&1
Write-Output '--- curl out ---'
Write-Output $out
