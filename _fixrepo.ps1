$t = (Get-Content -LiteralPath (Join-Path $env:TEMP 'ghflow_token.txt') -Raw).Trim()
$h = @{ Authorization = ("Bearer " + $t); Accept = 'application/vnd.github+json' }
$body = '{"name":"\u56fe\u8bba\u5b66\u4e60","description":"\u4ee3\u6570\u56fe\u8bba\u8bfe\u7a0b\u8d44\u6599\uff1a\u5168\u4e6613\u7ae0\u77e5\u8bc6\u70b9\u68b3\u7406 + \u7b2c8\u7ae0\u6388\u8bfe\u8bb2\u4e49\uff08MD/PDF\uff09","private":false,"auto_init":false}'
$bytes = [System.Text.Encoding]::UTF8.GetBytes($body)
$wantCps = '56fe8bba5b664e60'
$deadline = (Get-Date).AddMinutes(10)
$deleted = $false
while ((Get-Date) -lt $deadline) {
  if (-not $deleted) {
    try {
      Invoke-RestMethod -Uri 'https://api.github.com/repos/jhy28740-jhy/-' -Method Delete -Headers $h -TimeoutSec 25 | Out-Null
      $deleted = $true; Write-Output 'DELETED bogus'
    } catch {
      $resp = $_.Exception.Response
      if ($resp -and [int]$resp.StatusCode -eq 404) { $deleted = $true; Write-Output 'bogus already gone (404)' }
    }
  }
  if ($deleted) {
    try {
      $ok = Invoke-RestMethod -Uri 'https://api.github.com/user/repos' -Method Post -Headers $h -ContentType 'application/json; charset=utf-8' -Body $bytes -TimeoutSec 25
      if ($ok -and $ok.name) {
        $cps = ([int[]][char[]]$ok.name | ForEach-Object { $_.ToString('x4') }) -join ''
        Write-Output ("CREATED cps=" + $cps + " match=" + ($cps -eq $wantCps))
        Write-Output ("SSH=" + $ok.ssh_url)
        Write-Output ("HTML=" + $ok.html_url)
        exit 0
      }
    } catch {
      $resp = $_.Exception.Response
      if ($resp) {
        $sr = New-Object System.IO.StreamReader($resp.GetResponseStream()); $msg = $sr.ReadToEnd()
        if ($msg -match 'already exists') { Write-Output 'EXISTS-ALREADY'; exit 2 }
      }
    }
  }
  Start-Sleep -Seconds 8
}
Write-Output 'FIXREPO-FAIL timeout'
exit 1
