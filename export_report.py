"""
export_report.py
================
Převede REPORT.MD do samostatného REPORT.html (včetně GitHub-style alert boxů,
tabulek, KaTeX matematiky a vložených obrázků) a pomocí headless prohlížeče
(Chrome / Edge) vygeneruje tiskové REPORT.pdf ve vysokém rozlišení.
"""

import os
import re
import base64
import subprocess
from markdown_it import MarkdownIt

WORKSPACE = r"d:\prog\python\jupyter\Icaros48"
MD_PATH = os.path.join(WORKSPACE, "REPORT.MD")
HTML_PATH = os.path.join(WORKSPACE, "REPORT.html")
PDF_PATH = os.path.join(WORKSPACE, "REPORT.pdf")


def embed_images_as_base64(md_text: str, base_dir: str) -> str:
    """Nahradí relativní cesty k obrázkům v Markdownu za data:image/png;base64,..."""
    def repl(match):
        alt = match.group(1)
        rel_path = match.group(2)
        full_path = os.path.join(base_dir, rel_path.replace("/", os.sep))
        if os.path.exists(full_path):
            with open(full_path, "rb") as fp:
                b64 = base64.b64encode(fp.read()).decode("ascii")
            return f"![{alt}](data:image/png;base64,{b64})"
        return match.group(0)

    return re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", repl, md_text)


def convert_github_alerts(html_text: str) -> str:
    """Převede blockquote s [!IMPORTANT] / [!NOTE] na stylované GitHub alert boxy."""
    alert_Map = {
        "IMPORTANT": ("alert-important", "Důležité upozornění"),
        "NOTE": ("alert-note", "Poznámka"),
        "TIP": ("alert-tip", "Tip"),
        "WARNING": ("alert-warning", "Varování"),
        "CAUTION": ("alert-caution", "Pozor"),
    }
    for key, (css_cls, title) in alert_Map.items():
        pattern = r"<blockquote>\s*<p>\[!" + key + r"\]\s*"
        replacement = (
            f'<div class="gh-alert {css_cls}">'
            f'<div class="gh-alert-title">{title}</div><p>'
        )
        # Nahradíme otevírací blockquote a příslušný uzavírací </blockquote>
        while re.search(pattern, html_text):
            m = re.search(pattern, html_text)
            start = m.start()
            end_bq = html_text.find("</blockquote>", m.end())
            if end_bq != -1:
                inner = html_text[m.end():end_bq]
                html_text = (
                    html_text[:start]
                    + replacement
                    + inner
                    + "</div>"
                    + html_text[end_bq + len("</blockquote>"):]
                )
            else:
                break
    return html_text


def protect_math(md_text: str):
    """Dočasně ochrání LaTeX výrazy $...$ a $$...$$ před markdown parserem."""
    placeholders = {}
    counter = [0]

    def store(m):
        key = f"@@MATHBLOCK{counter[0]}@@"
        counter[0] += 1
        placeholders[key] = m.group(0)
        return key

    # Display math $$...$$
    md_text = re.sub(r"\$\$[\s\S]+?\$\$", store, md_text)
    # Inline math $...$
    md_text = re.sub(r"(?<!\\)\$[^\$\n]+?\$", store, md_text)
    return md_text, placeholders


def restore_math(html_text: str, placeholders: dict) -> str:
    for k, v in placeholders.items():
        html_text = html_text.replace(k, v)
    return html_text


def build_html() -> None:
    with open(MD_PATH, "r", encoding="utf-8") as fp:
        md_raw = fp.read()

    md_embedded = embed_images_as_base64(md_raw, WORKSPACE)
    md_protected, math_map = protect_math(md_embedded)

    md = MarkdownIt("commonmark", {"html": True, "typographer": False}).enable("table")
    body_html = md.render(md_protected)
    body_html = restore_math(body_html, math_map)
    body_html = convert_github_alerts(body_html)

    full_html = f"""<!DOCTYPE html>
<html lang="cs">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>REPORT - Analýza měřených dat a dokumentace programu</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/contrib/auto-render.min.js"
    onload="renderMathInElement(document.body, {{
        delimiters: [
            {{left: '$$', right: '$$', display: true}},
            {{left: '$', right: '$', display: false}}
        ],
        throwOnError: false
    }});"></script>
<style>
    @page {{
        size: A4 portrait;
        margin: 14mm 14mm 16mm 14mm;
    }}
    * {{
        box-sizing: border-box;
        -webkit-print-color-adjust: exact !important;
        print-color-adjust: exact !important;
    }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        font-size: 10.5pt;
        line-height: 1.55;
        color: #1f2328;
        background: #ffffff;
        max-width: 980px;
        margin: 0 auto;
        padding: 20px 28px;
    }}
    @media print {{
        body {{
            max-width: 100%;
            padding: 0;
            font-size: 9.5pt;
        }}
        h2, h3 {{
            page-break-after: avoid;
            break-after: avoid;
        }}
        table, figure, img, .gh-alert, pre {{
            page-break-inside: avoid;
            break-inside: avoid;
        }}
        tr {{
            page-break-inside: avoid;
            break-inside: avoid;
        }}
    }}
    h1 {{
        font-size: 19pt;
        border-bottom: 2px solid #d0d7de;
        padding-bottom: 8px;
        margin-top: 0;
        color: #0d47a1;
    }}
    h2 {{
        font-size: 13.5pt;
        border-bottom: 1px solid #d8dee4;
        padding-bottom: 5px;
        margin-top: 22px;
        color: #1565c0;
    }}
    h3 {{
        font-size: 11.5pt;
        margin-top: 16px;
        margin-bottom: 6px;
        color: #263238;
    }}
    hr {{
        border: 0;
        border-top: 1px solid #d0d7de;
        margin: 18px 0;
    }}
    a {{
        color: #0969da;
        text-decoration: none;
    }}
    code {{
        font-family: Consolas, "Courier New", monospace;
        font-size: 0.92em;
        background-color: #f6f8fa;
        padding: 2px 5px;
        border-radius: 4px;
        border: 1px solid #eaeef2;
    }}
    pre {{
        background-color: #f6f8fa;
        border: 1px solid #d0d7de;
        border-radius: 6px;
        padding: 12px 14px;
        overflow-x: auto;
    }}
    pre code {{
        background: none;
        border: none;
        padding: 0;
        font-size: 9pt;
    }}
    table {{
        width: 100%;
        border-collapse: collapse;
        margin: 12px 0 16px 0;
        font-size: 8.8pt;
    }}
    th, td {{
        border: 1px solid #d0d7de;
        padding: 5px 8px;
        text-align: left;
        vertical-align: top;
    }}
    th {{
        background-color: #f0f4f8;
        font-weight: 600;
        color: #0d47a1;
    }}
    tr:nth-child(even) td {{
        background-color: #fafbfc;
    }}
    img {{
        display: block;
        max-width: 100%;
        height: auto;
        margin: 10px auto 16px auto;
        border: 1px solid #d0d7de;
        border-radius: 6px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    }}
    .gh-alert {{
        padding: 10px 14px;
        margin: 12px 0;
        border-left: 4px solid #0969da;
        border-radius: 4px;
        background-color: #f6f8fa;
    }}
    .gh-alert p, .gh-alert ul {{
        margin: 4px 0;
    }}
    .gh-alert-title {{
        font-weight: 700;
        margin-bottom: 4px;
        font-size: 10pt;
    }}
    .alert-important {{
        border-left-color: #8250df;
        background-color: #fbfaff;
    }}
    .alert-important .gh-alert-title {{
        color: #8250df;
    }}
    .alert-note {{
        border-left-color: #0969da;
        background-color: #f4f9ff;
    }}
    .alert-note .gh-alert-title {{
        color: #0969da;
    }}
</style>
</head>
<body>
{body_html}
</body>
</html>
"""
    with open(HTML_PATH, "w", encoding="utf-8") as fp:
        fp.write(full_html)
    print(f"HTML úspěšně vygenerováno: {HTML_PATH} ({os.path.getsize(HTML_PATH) / 1024:.1f} KB)")


def build_pdf() -> None:
    browsers = [
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    ]
    browser_exe = next((b for b in browsers if os.path.exists(b)), None)
    if not browser_exe:
        raise RuntimeError("Nebyl nalezen prohlížeč Chrome ani Edge pro generování PDF.")

    cmd = [
        browser_exe,
        "--headless",
        "--disable-gpu",
        "--no-pdf-header-footer",
        "--print-to-pdf-no-header",
        "--virtual-time-budget=5000",
        f"--print-to-pdf={PDF_PATH}",
        f"file:///{HTML_PATH.replace(os.sep, '/')}",
    ]
    subprocess.run(cmd, check=True)
    print(f"PDF úspěšně vygenerováno: {PDF_PATH} ({os.path.getsize(PDF_PATH) / 1024:.1f} KB)")


if __name__ == "__main__":
    build_html()
    build_pdf()
