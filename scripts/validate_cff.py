#!/usr/bin/env python3
# scripts/validate_cff.py
import sys
import yaml
from pathlib import Path

def validate_and_print(path):
    p = Path(path)
    if not p.exists():
        print("File not found:", path); return 2
    data = yaml.safe_load(p.read_text())
    # Basic checks
    title = data.get("title") or data.get("name") or "<no title>"
    authors = data.get("authors") or []
    date = data.get("date-released") or data.get("date") or "<no date>"
    version = data.get("version", "<no version>")
    print("Title :", title)
    print("Version:", version)
    print("Date  :", date)
    print("Authors:")
    for a in authors:
        fn = a.get("family-names") or a.get("family") or ""
        gn = a.get("given-names") or a.get("given") or ""
        email = a.get("email") or ""
        aff = a.get("affiliation") or ""
        print(f" - {gn} {fn}".strip() + (f" <{email}>" if email else "") + (f" [{aff}]" if aff else ""))
    # done
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scripts/validate_cff.py CITATION.cff")
        sys.exit(2)
    sys.exit(validate_and_print(sys.argv[1]))
