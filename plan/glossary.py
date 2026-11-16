#!/usr/bin/env python3
"""Rebuild the Parts 2+ section of dictionary.qmd from plan/glossary/<slug>.md.

Part 1 is hand-written and kept. Everything after the marker
'<!-- BEGIN GENERATED -->' is regenerated. Run: python3 plan/glossary.py
"""
import os, re
rows = []
for line in open('plan/units.md'):
    m = re.match(r'\| (\d+) \| `([a-z0-9-]+)` \| (.*?) \| (Part \d+): (.*?) \|', line)
    if m: rows.append((int(m.group(1)), m.group(2), m.group(3), m.group(4), m.group(5)))
out, cur = [], None
for n, slug, title, part, pname in sorted(rows):
    if n <= 4: continue
    f = f'plan/glossary/{slug}.md'
    if not os.path.exists(f): continue
    terms = []
    for l in open(f):
        if l.startswith('## Install'): break
        l = l.strip()
        if l.startswith('- **'): l = l[2:]
        if l.startswith('**'): terms.append(l)
    if not terms: continue
    if part != cur:
        cur = part
        out += ['', f'## {part}: {pname}']
    out += ['', f'### {title}', '', f'See [{title}]({slug}/lec.qmd).', '']
    out += [t + '\n' for t in terms]
s = open('dictionary.qmd').read()
mark = '<!-- BEGIN GENERATED (plan/glossary.py) -->'
if mark not in s:
    i = s.index('<!-- Terms to add')
    s = s[:i] + mark + '\n'
else:
    s = s[:s.index(mark) + len(mark)] + '\n'
s += '\n'.join(out).replace('\n\n\n', '\n\n') + '\n'
open('dictionary.qmd', 'w').write(s)
print(sum(1 for l in out if l.startswith('**')), 'terms')
