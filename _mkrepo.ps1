$t = (Get-Content -LiteralPath (Join-Path $env:TEMP 'ghflow_token.txt') -Raw).Trim()
$body = '{"name":"图论学习","description":"代数图论课程资料：全书13章知识点梳理 + 第8章授课讲义（MD/PDF）","private":false,"auto_init":false}'
$bytes = [System.Text.Encoding]::UTF8.GetBytes($body)
$h = @{ Authorization = ("Bearer " + $t); Accept = 'application/vnd.github+json' }
$ok = $null
for ($i = 1; $i -le 6; $i++) {
  try {
    $ok = Invoke-RestMethod -Uri 'https://api.github.com/user/repos' -Method Post -Headers $h -ContentType 'application/json; charset=utf-8' -Body $bytes -TimeoutSec 30
    break
  } catch {
    $resp = $_.Exception.Response
    if ($resp) {
      $sr = New-Object System.IO.StreamReader($resp.GetResponseStream())
      $msg = $sr.ReadToEnd()
      Write-Output ("HTTP " + [int]$resp.StatusCode + ": " + $msg)
      if ($msg -notmatch 'Problems parsing JSON') { break }
    } else { Write-Output ("NET-ERR try " + $i + ": " + $_.Exception.Message) }
    Start-Sleep -Seconds 3
  }
}
if ($ok -and $ok.full_name) {
  Write-Output ("REPO-CREATED=" + $ok.full_name)
  Write-Output ("HTML=" + $ok.html_url)
  Write-Output ("SSH=" + $ok.ssh_url)
} else { Write-Output 'REPO-FAIL' }
