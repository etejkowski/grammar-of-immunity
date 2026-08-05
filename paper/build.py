#!/usr/bin/env python3
"""
Build submission-ready manuscript files
=======================================

Produces, in paper/build/:
    manuscript-submission.md    manuscript with figures embedded inline
    manuscript.docx             for journal portals (Editorial Manager et al.)
    manuscript.html             for quick review in a browser
    manuscript.pdf              for OSF Preprints, which wants a PDF

The source manuscript names its figures in the captions but does not embed
them, which keeps the Markdown readable on GitHub. This script inserts the
image references at each caption before handing the file to pandoc, so the
source stays clean and the deliverable is complete.

Requires pandoc (brew install pandoc). No LaTeX needed: the PDF is printed
from the HTML by headless Chrome, which is already on this machine. If Chrome
is absent the PDF step is skipped with a warning and the other files still
build.

Usage:
    python3 paper/build.py
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'manuscript.md')
BUILD = os.path.join(HERE, 'build')
FIGDIR = os.path.join(HERE, 'figures')

CAPTION = re.compile(r'^\*\*Figure (\d+)\*\* \(`([^`]+)`\)\.', re.M)

CHROME = ('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')

# Print stylesheet for the PDF. Letter portrait with 1 inch margins to match
# the DOCX, figures and tables kept off page breaks, and tables set small
# enough that the widest one (Table 5) fits the text block.
PRINT_CSS = """
@page { size: letter portrait; margin: 1in; }
html { font-size: 11pt; }
body { font-family: Georgia, "Times New Roman", serif; line-height: 1.45;
       max-width: none; margin: 0; padding: 0; color: #000; }
h1 { font-size: 16pt; line-height: 1.25; }
h2 { font-size: 13pt; margin-top: 1.4em; }
h3 { font-size: 11.5pt; }
h1, h2, h3 { break-after: avoid; page-break-after: avoid; }
p { orphans: 3; widows: 3; }
img { max-width: 100%; height: auto; }
figure, table { break-inside: avoid; page-break-inside: avoid; }
table { border-collapse: collapse; font-size: 8.5pt; width: 100%;
        margin: 0.8em 0; }
th, td { border: 1px solid #999; padding: 2px 5px; text-align: left; }
th { background: #f0f0f0; }
code, pre { font-family: Menlo, Consolas, monospace; font-size: 8.5pt; }
a { color: #000; text-decoration: none; }
hr { display: none; }
"""


def make_pdf(html_path, pdf_path):
    """
    Print the HTML build to PDF with headless Chrome.

    OSF Preprints, Preprints.org and Zenodo all want a PDF, and pandoc's own
    PDF writers need a LaTeX or WebKit engine that is not installed here.
    Chrome renders the same HTML the browser shows, so what you proof is what
    gets deposited. Returns True on success.
    """
    if not os.path.exists(CHROME):
        print("  Chrome not found; skipping PDF. Install Chrome, or export the")
        print("  DOCX to PDF from Pages, to produce paper/build/manuscript.pdf")
        return False
    if os.path.exists(pdf_path):
        os.remove(pdf_path)
    cmd = [CHROME, '--headless=new', '--disable-gpu', '--no-pdf-header-footer',
           '--run-all-compositor-stages-before-draw',
           '--virtual-time-budget=30000',
           f'--print-to-pdf={pdf_path}', 'file://' + html_path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(pdf_path):
        print(f"FAILED manuscript.pdf:\n{r.stderr[-600:]}")
        return False
    with open(pdf_path, 'rb') as f:
        pages = len(re.findall(rb'/Type\s*/Page[^s]', f.read()))
    print(f"  wrote {os.path.relpath(pdf_path)}  "
          f"({os.path.getsize(pdf_path):,} bytes, {pages} pages)")
    return True


def set_letter_page_size(path):
    """
    Force 8.5 x 11 inch portrait with 1 inch margins in a built .docx.

    bioRxiv and most publishers require US Letter portrait for reliable PDF
    conversion, but pandoc does not propagate the reference document's page
    size into its output, so the file would otherwise inherit whatever default
    the reading application picks. Sizes are in twips: 12240 x 15840 = 8.5 x 11.
    """
    import zipfile
    name = 'word/document.xml'
    with zipfile.ZipFile(path) as z:
        items = {n: z.read(n) for n in z.namelist()}
    xml = items[name].decode('utf-8')
    if '<w:pgSz' in xml:
        return False
    pg = ('<w:pgSz w:w="12240" w:h="15840" w:orient="portrait"/>'
          '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
          'w:header="720" w:footer="720" w:gutter="0"/>')
    if re.search(r'<w:sectPr\b[^>]*>', xml):
        xml = re.sub(r'(<w:sectPr\b[^>]*>)', r'\1' + pg, xml, count=1)
    else:
        xml = xml.replace('</w:body>', f'<w:sectPr>{pg}</w:sectPr></w:body>')
    items[name] = xml.encode('utf-8')
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for n, data in items.items():
            z.writestr(n, data)
    return True


def write_plain_abstract(text, dest):
    """
    Emit the Abstract as plain text for preprint and journal submission forms.

    Every portal wants the abstract pasted as text, and none of them accept
    Markdown emphasis. Structured headings ("Background.", "Methods.") are kept
    because they are part of the abstract's content; only the bold markers go.
    Returns the word count, which portals cap and which the venue files quote.
    """
    m = re.search(r'\n## Abstract\n(.*?)\n---\n', text, re.S)
    if not m:
        print("  WARNING: no Abstract section found; skipping plain abstract")
        return 0
    body = m.group(1).strip()
    body = re.sub(r'\*\*(.+?)\*\*', r'\1', body)     # drop bold markers
    body = re.sub(r'\[[\d,\u2013\u2014 -]+\]', '', body)  # drop citation brackets
    body = re.sub(r' +([.,;])', r'\1', body)
    body = re.sub(r'[ \t]+\n', '\n', body)
    with open(dest, 'w', encoding='utf-8') as f:
        f.write(body + '\n')
    words = len([w for w in re.split(r'\s+', re.sub(
        r'(?s)\nKeywords:.*', '', body)) if w])
    print(f"  wrote {os.path.relpath(dest)}  ({words} words, excluding keywords)")
    return words


def main():
    if not shutil.which('pandoc'):
        print("pandoc not found. Install with: brew install pandoc")
        return 1
    os.makedirs(BUILD, exist_ok=True)

    text = open(SRC, encoding='utf-8').read()

    # Insert the image immediately before each caption paragraph.
    missing = []
    def repl(m):
        num, fname = m.group(1), m.group(2)
        if not os.path.exists(os.path.join(FIGDIR, fname)):
            missing.append(fname)
            return m.group(0)
        return f'![](figures/{fname}){{width=100%}}\n\n' + m.group(0)

    text, n = CAPTION.subn(repl, text)
    print(f"embedded {n} figures")
    if missing:
        print("MISSING figure files:", ', '.join(missing))
        return 1

    out_md = os.path.join(BUILD, 'manuscript-submission.md')
    with open(out_md, 'w', encoding='utf-8') as f:
        f.write(text)

    # pandoc runs from paper/ so the relative figures/ paths resolve
    common = ['--from', 'markdown+pipe_tables+tex_math_dollars',
              '--resource-path', HERE, '--standalone']
    # reference.docx is pandoc's default template with the Microsoft fonts
    # (Aptos, Consolas, Times New Roman, Segoe UI) swapped for macOS-native
    # equivalents, so Pages and Preview stop reporting missing fonts and the
    # PDF export is not silently substituted. Regenerate with:
    #   pandoc --print-default-data-file reference.docx > ref.docx
    # then patch the font names in word/styles.xml and word/theme/theme1.xml.
    ref = os.path.join(HERE, 'reference.docx')
    docx_extra = ['--reference-doc', ref] if os.path.exists(ref) else []

    # The print HTML is the PDF's source: same content, no table of contents,
    # and the print stylesheet above instead of pandoc's screen defaults.
    css = os.path.join(BUILD, 'print.css')
    with open(css, 'w', encoding='utf-8') as f:
        f.write(PRINT_CSS)

    targets = [('manuscript.docx', docx_extra),
               ('manuscript.html', ['--embed-resources', '--toc']),
               ('manuscript-print.html', ['--embed-resources', '--css', css])]
    for name, extra in targets:
        dest = os.path.join(BUILD, name)
        cmd = ['pandoc', out_md, '-o', dest] + common + extra
        r = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAILED {name}:\n{r.stderr[:600]}")
            return 1
        if name.endswith('.docx') and set_letter_page_size(dest):
            print("    set page size to 8.5 x 11 in portrait")
        size = os.path.getsize(dest)
        print(f"  wrote {os.path.relpath(dest)}  ({size:,} bytes)")

    make_pdf(os.path.join(BUILD, 'manuscript-print.html'),
             os.path.join(BUILD, 'manuscript.pdf'))
    write_plain_abstract(text, os.path.join(BUILD, 'abstract-plain.txt'))

    print("\nFigures are embedded. Most journals also want figures as separate")
    print("files at submission; paper/figures/*.png serves that directly.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
