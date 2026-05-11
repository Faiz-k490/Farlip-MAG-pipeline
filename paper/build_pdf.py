"""Render paper/paper.md to paper/paper.pdf via Chrome headless.

Usage:
    .venv/bin/python paper/build_pdf.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent
MD_PATH = ROOT / "paper.md"
HTML_PATH = ROOT / "paper.html"
PDF_PATH = ROOT / "paper.pdf"

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

CSS = """
@page {
  size: Letter;
  margin: 1.0in 1.1in 1.0in 1.1in;
  @bottom-center {
    content: counter(page);
    font-family: "Times New Roman", Times, serif;
    font-size: 9.5pt;
    color: #333;
  }
}
html { font-size: 10.5pt; }
body {
  font-family: "Times New Roman", "Times", "Liberation Serif", serif;
  line-height: 1.42;
  color: #000;
  margin: 0;
}
h1, h2, h3 {
  font-family: "Times New Roman", "Times", serif;
  font-weight: bold;
  page-break-after: avoid;
}
h1 {
  font-size: 1.18rem;
  margin: 1.5em 0 0.5em 0;
}
h2 {
  font-size: 1.06rem;
  margin: 1.1em 0 0.35em 0;
}
h3 {
  font-size: 1.0rem;
  font-style: italic;
  font-weight: normal;
  margin: 0.9em 0 0.3em 0;
}
p {
  margin: 0 0 0.55em 0;
  text-align: justify;
  hyphens: auto;
  text-indent: 0;
}
/* Indent paragraphs after the first within a section, classical style */
h1 + p, h2 + p, h3 + p, hr + p { text-indent: 0; }
p + p { text-indent: 1.5em; }
ul, ol { margin: 0.3em 0 0.6em 1.3em; }
li { margin: 0.1em 0; }
code {
  font-family: "Courier New", "Courier", monospace;
  font-size: 0.9em;
}
pre {
  font-family: "Courier New", "Courier", monospace;
  font-size: 0.84em;
  line-height: 1.35;
  margin: 0.7em 1.2em;
  white-space: pre-wrap;
  page-break-inside: avoid;
}
pre code { font-size: inherit; }
table {
  border-collapse: collapse;
  margin: 0.7em auto;
  font-size: 0.93em;
  page-break-inside: avoid;
}
th, td {
  border-top: 0.5pt solid #000;
  border-bottom: 0.5pt solid #000;
  padding: 3px 10px;
  text-align: left;
  vertical-align: top;
}
th { font-weight: bold; border-bottom: 1pt solid #000; }
hr { border: none; margin: 1.0em 0; }
.frontmatter { margin-bottom: 1.6em; text-align: center; }
.frontmatter .title {
  display: block;
  font-size: 1.35rem;
  font-weight: bold;
  margin-bottom: 0.6em;
  line-height: 1.3;
}
.frontmatter .meta { font-size: 0.95em; line-height: 1.5; }
"""

HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
{frontmatter_html}
{body}
</body>
</html>
"""


def parse_frontmatter(text: str) -> tuple[dict, str]:
    """Pull out YAML-ish frontmatter delimited by --- ... --- at the top."""
    if not text.startswith("---\n"):
        return {}, text
    end = text.find("\n---\n", 4)
    if end == -1:
        return {}, text
    raw = text[4:end]
    body = text[end + 5 :]
    meta: dict = {}
    current_key = None
    for line in raw.splitlines():
        if not line.strip():
            continue
        if line.startswith("  - ") and current_key:
            meta.setdefault(current_key, []).append(line[4:].strip())
        elif ":" in line and not line.startswith(" "):
            key, _, val = line.partition(":")
            key = key.strip()
            val = val.strip().strip('"')
            current_key = key
            if val:
                meta[key] = val
            else:
                meta[key] = []
    return meta, body


def render_frontmatter(meta: dict) -> str:
    if not meta:
        return ""
    parts = []
    title = meta.get("title", "")
    if title:
        parts.append(f'<span class="title">{title}</span>')
    meta_lines = []
    authors = meta.get("authors")
    if isinstance(authors, list) and authors:
        meta_lines.append(", ".join(authors))
    if meta.get("affiliation"):
        meta_lines.append(f'<i>{meta["affiliation"]}</i>')
    if meta.get("date"):
        meta_lines.append(meta["date"])
    if meta_lines:
        parts.append('<div class="meta">' + "<br>".join(meta_lines) + "</div>")
    return f'<div class="frontmatter">{"".join(parts)}</div>'


def main() -> int:
    if not MD_PATH.exists():
        print(f"missing: {MD_PATH}", file=sys.stderr)
        return 2
    if not Path(CHROME).exists():
        print(f"missing Chrome at {CHROME}", file=sys.stderr)
        return 2

    text = MD_PATH.read_text(encoding="utf-8")
    meta, body_md = parse_frontmatter(text)

    md = markdown.Markdown(
        extensions=["tables", "fenced_code", "codehilite", "toc", "sane_lists"],
        extension_configs={"codehilite": {"guess_lang": False, "noclasses": True}},
    )
    body_html = md.convert(body_md)

    html = HTML_TEMPLATE.format(
        title=meta.get("title", "paper"),
        css=CSS,
        frontmatter_html=render_frontmatter(meta),
        body=body_html,
    )
    HTML_PATH.write_text(html, encoding="utf-8")
    print(f"wrote {HTML_PATH}")

    cmd = [
        CHROME,
        "--headless=new",
        "--disable-gpu",
        "--no-pdf-header-footer",
        "--no-margins",
        "--virtual-time-budget=10000",
        f"--print-to-pdf={PDF_PATH}",
        HTML_PATH.as_uri(),
    ]
    print("running:", " ".join(shutil.which("echo") and a or a for a in cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("chrome stdout:", result.stdout, file=sys.stderr)
        print("chrome stderr:", result.stderr, file=sys.stderr)
        return result.returncode
    if not PDF_PATH.exists():
        print("PDF not produced", file=sys.stderr)
        return 1
    size_kb = PDF_PATH.stat().st_size / 1024
    print(f"wrote {PDF_PATH} ({size_kb:.1f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
