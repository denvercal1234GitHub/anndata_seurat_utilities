#!/usr/bin/env python3
# scripts/cff_to_bib.py
import sys, yaml, re
from pathlib import Path
from datetime import datetime

def make_bibkey(authors, year, title):
    # simple key: FirstAuthorFamilyYYYYshorttitle
    if authors:
        fam = authors[0].get("family-names") or authors[0].get("family") or "anon"
    else:
        fam = "anon"
    fam = re.sub(r'[^A-Za-z0-9]', '', fam)
    short = re.sub(r'[^A-Za-z0-9]', '', (title.split()[0] if title else "pkg"))
    return f"{fam}{year}{short}"

def format_authors_bibtex(authors):
    if not authors:
        return ""
    parts = []
    for a in authors:
        gn = a.get("given-names") or a.get("given") or ""
        fn = a.get("family-names") or a.get("family") or ""
        if gn and fn:
            parts.append(f"{gn} {fn}")
        elif fn:
            parts.append(fn)
        elif gn:
            parts.append(gn)
    return " and ".join(parts)

def to_bib(path):
    p = Path(path)
    data = yaml.safe_load(p.read_text())

    title = data.get("title", "")
    authors = data.get("authors", [])
    year = None
    if data.get("date-released"):
        year = data["date-released"][:4]
    elif data.get("date"):
        year = data["date"][:4]
    else:
        year = str(datetime.now().year)
    version = data.get("version", "")
    url = data.get("url") or data.get("repository-code") or ""
    note = data.get("note", "")
    bibkey = make_bibkey(authors, year, title)
    author_field = format_authors_bibtex(authors)

    bib = []
    bib.append(f"@software{{{bibkey},")
    if author_field:
        bib.append(f"  author = {{{author_field}}},")
    bib.append(f"  title = {{{title}}},")
    bib.append(f"  year = {{{year}}},")
    if url:
        bib.append(f"  url = {{{url}}},")
    if version:
        bib.append(f"  version = {{{version}}},")
    if note:
        bib.append(f"  note = {{{note}}},")
    bib.append("}")
    return "\n".join(bib)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scripts/cff_to_bib.py CITATION.cff [output.bib]")
        sys.exit(2)
    out = None
    if len(sys.argv) > 2:
        out = Path(sys.argv[2])
    bib = to_bib(sys.argv[1])
    if out:
        out.write_text(bib)
        print("Wrote", out)
    else:
        print(bib)
