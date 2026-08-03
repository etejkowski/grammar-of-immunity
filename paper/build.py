#!/usr/bin/env python3
"""
Build submission-ready manuscript files
=======================================

Produces, in paper/build/:
    manuscript-submission.md    manuscript with figures embedded inline
    manuscript.docx             for bioRxiv and most journal portals
    manuscript.html             for quick review in a browser

The source manuscript names its figures in the captions but does not embed
them, which keeps the Markdown readable on GitHub. This script inserts the
image references at each caption before handing the file to pandoc, so the
source stays clean and the deliverable is complete.

Requires pandoc (brew install pandoc). No LaTeX needed for DOCX or HTML.

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
    targets = [('manuscript.docx', docx_extra),
               ('manuscript.html', ['--embed-resources', '--toc'])]
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

    print("\nFigures are embedded. Most journals also want figures as separate")
    print("files at submission; paper/figures/*.png serves that directly.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
