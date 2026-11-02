import http from 'node:http';
import fs from 'node:fs';
import path from 'node:path';
import net from 'node:net';
import { pathToFileURL } from 'node:url';

const ROOT = 'C:/Users/hyji11/.agents/skills/dashi-ppt-skill/skills/dashi-ppt/project';
const deckDir = path.resolve(process.argv[2]);
const outFile = path.resolve(process.argv[3]);
const title = process.argv[4] || 'Deck PDF Export';

const MIME = {
  '.html': 'text/html;charset=utf-8', '.js': 'text/javascript;charset=utf-8', '.mjs': 'text/javascript;charset=utf-8',
  '.css': 'text/css;charset=utf-8', '.json': 'application/json;charset=utf-8', '.map': 'application/json;charset=utf-8',
  '.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg', '.svg': 'image/svg+xml', '.ico': 'image/x-icon',
  '.woff': 'font/woff', '.woff2': 'font/woff2', '.ttf': 'font/ttf', '.otf': 'font/otf',
  '.mp4': 'video/mp4', '.webm': 'video/webm', '.avif': 'image/avif', '.txt': 'text/plain;charset=utf-8'
};

const server = http.createServer((req, res) => {
  try {
    const u = new URL(req.url, 'http://127.0.0.1');
    let p = decodeURIComponent(u.pathname);
    if (p === '/') p = '/index.html';
    const fp = path.join(deckDir, p);
    if (!fp.startsWith(deckDir) || !fs.existsSync(fp) || fs.statSync(fp).isDirectory()) {
      res.writeHead(404, { 'content-type': 'text/plain;charset=utf-8' }); res.end('not found'); return;
    }
    res.writeHead(200, { 'content-type': MIME[path.extname(fp).toLowerCase()] || 'application/octet-stream', 'cache-control': 'no-store' });
    fs.createReadStream(fp).pipe(res);
  } catch (e) { res.writeHead(500); res.end(String(e)); }
});

const port = await new Promise((resolve, reject) => {
  const probe = net.createServer();
  probe.on('error', reject);
  probe.listen(0, '127.0.0.1', () => { const q = probe.address().port; probe.close(() => resolve(q)); });
});
await new Promise(r => server.listen(port, '127.0.0.1', r));
console.log('static preview on http://127.0.0.1:' + port + '/');

const pw = await import(pathToFileURL(ROOT + '/node_modules/playwright-core/index.js').href);
const chromium = pw.chromium || (pw.default && pw.default.chromium);
const { launchExportBrowser } = await import(pathToFileURL(ROOT + '/scripts/preview/launch-export-browser.mjs').href);
const { exportScreenshotPdfFromUrl } = await import(pathToFileURL(ROOT + '/packages/html-deck-to-pptx/dist/screenshot.min.mjs').href);

let browser = null;
try {
  browser = await launchExportBrowser(chromium, { fallbackTmpDirs: [path.join(path.dirname(outFile), '.browser-tmp')], log: m => console.warn(m) });
  const result = await exportScreenshotPdfFromUrl(browser, `http://127.0.0.1:${port}/`, { outFile, title });
  console.log('PDF exported: ' + outFile + ' (' + result.pages + ' page(s))');
} finally {
  if (browser) await browser.close().catch(() => {});
  server.close();
}

