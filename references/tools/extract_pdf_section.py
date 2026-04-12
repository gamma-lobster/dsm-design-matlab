#!/usr/bin/env python3
"""Search a PDF and print pages matching a text pattern.

Example:
    python3 references/extract_pdf_section.py \
        "UNDERSTANDING DELTA-SIGMA DATA CONVERTERS 2nd Edition.pdf" \
        --pattern "Excess Loop Delay" --context-pages 2
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import fitz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", help="Path to the PDF file")
    parser.add_argument(
        "--pattern",
        required=True,
        help="Regex pattern to search for in extracted text",
    )
    parser.add_argument(
        "--context-pages",
        type=int,
        default=0,
        help="Include this many pages before and after each matching page",
    )
    parser.add_argument(
        "--max-pages",
        type=int,
        default=20,
        help="Maximum number of pages to print",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pdf_path = Path(args.pdf)
    if not pdf_path.exists():
        print(f"PDF not found: {pdf_path}", file=sys.stderr)
        return 1

    pattern = re.compile(args.pattern, re.IGNORECASE)
    doc = fitz.open(pdf_path)

    matches = []
    for i in range(doc.page_count):
        text = doc.load_page(i).get_text("text")
        if pattern.search(text):
            start = max(0, i - args.context_pages)
            end = min(doc.page_count - 1, i + args.context_pages)
            matches.extend(range(start, end + 1))

    pages = []
    seen = set()
    for idx in matches:
        if idx not in seen:
            seen.add(idx)
            pages.append(idx)
        if len(pages) >= args.max_pages:
            break

    if not pages:
        print("No matches found.", file=sys.stderr)
        return 2

    print(f"PDF: {pdf_path}")
    print(f"Pages: {doc.page_count}")
    print(f"Pattern: {args.pattern}")
    print()

    for idx in pages:
        text = doc.load_page(idx).get_text("text").strip()
        print(f"===== PAGE {idx + 1} =====")
        print(text)
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
