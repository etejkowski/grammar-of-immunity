#!/usr/bin/env python3
"""
Check paper/highlights.md against Elsevier's limit
==================================================

Elsevier allows 3-5 highlights of at most 85 characters each, including spaces.
Exceeding the limit is a common desk-return reason, so this verifies it rather
than trusting a manual count.

Usage:
    python3 paper/check_highlights.py
"""

import os
import re
import sys

LIMIT = 85
MAX_BULLETS = 5
MIN_BULLETS = 3
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'highlights.md')


def main():
    text = open(SRC, encoding='utf-8').read()
    # the bullet list, not the character-count table
    bullets = [m.group(1).strip() for m in
               re.finditer(r'^- (.+)$', text, re.M)
               if not m.group(1).startswith(('Bullet', 'If the journal',
                                             'Avoid the'))]
    if not bullets:
        print('no bullets found')
        return 1

    ok = True
    print(f"{'#':>2} {'chars':>5}  bullet")
    for i, b in enumerate(bullets, 1):
        n = len(b)
        flag = '' if n <= LIMIT else f'  <-- OVER by {n - LIMIT}'
        if n > LIMIT:
            ok = False
        print(f"{i:>2} {n:>5}  {b}{flag}")

    if not MIN_BULLETS <= len(bullets) <= MAX_BULLETS:
        print(f"\ncount {len(bullets)} outside Elsevier's "
              f"{MIN_BULLETS}-{MAX_BULLETS}")
        ok = False

    # the table in highlights.md must agree with the real counts
    table = {int(m.group(1)): int(m.group(2)) for m in
             re.finditer(r'^\| (\d) \| (\d+) \|', text, re.M)}
    for i, b in enumerate(bullets, 1):
        if i in table and table[i] != len(b):
            print(f"\ntable says bullet {i} is {table[i]} chars, "
                  f"actual {len(b)}")
            ok = False

    print('\nPASS' if ok else '\nFAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
