$p = "第8章课程讲义_强正则图与设计.md"
$bytes = [System.IO.File]::ReadAllBytes($p)
$hasBom = ($bytes.Length -ge 3 -and $bytes[0] -eq 0xEF -and $bytes[1] -eq 0xBB -and $bytes[2] -eq 0xBF)
$t = [System.IO.File]::ReadAllText($p, [System.Text.UTF8Encoding]::new($false))
$pairs = @(
  @('=\#\{\text{$u$ 与 $(1,1)$ 的公共邻居}\}=a=n .', '=\#\{u\ \text{与}\ (1,1)\ \text{的公共邻居}\}=a=n .'),
  @('\Longleftrightarrow\ \text{$s$ 是 2 阶元}.',      '\Longleftrightarrow\ s\ \text{是 2 阶元}.'),
  @('\cdot(\text{$p_i$ 的首项系数})\ne0',             '\cdot(p_i\ \text{的首项系数})\ne0')
)
foreach ($pr in $pairs) {
  if ($t.Contains($pr[0])) { $t = $t.Replace($pr[0], $pr[1]); Write-Output "OK  : $($pr[0].Substring(0,[Math]::Min(46,$pr[0].Length)))" }
  else { Write-Output "MISS: $($pr[0].Substring(0,[Math]::Min(46,$pr[0].Length)))" }
}
$enc = [System.Text.UTF8Encoding]::new($hasBom)
[System.IO.File]::WriteAllText($p, $t, $enc)
Write-Output "bom=$hasBom written"
