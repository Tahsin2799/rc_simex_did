#!/usr/bin/env bash
# Build the arXiv source bundle from the rendered manuscript and test that it
# compiles on its own.
#   bash paper1_rc_simex_did/submission/make_arxiv_bundle.sh     (after rendering paper_all_in_one.Rmd)
# Needs pdflatex and bibtex on PATH (TinyTeX). Temporary build directories go
# under $TMPDIR. Writes arxiv/ and arxiv_bundle.tar.gz next to this script.
# arXiv does not run BibTeX, so the bundle carries ms.bbl and no .bib file.
set -euo pipefail
cd "$(dirname "$0")"
OUT=arxiv
rm -rf "$OUT" arxiv_bundle.tar.gz
mkdir -p "$OUT/figures"
cp paper_all_in_one.tex "$OUT/ms.tex"
for f in $(grep -o 'includegraphics[^{]*{[^}]*}' paper_all_in_one.tex | sed 's/.*{\(.*\)}/\1/'); do
  cp "$f.pdf" "$OUT/$f.pdf"
done

latex3() { pdflatex -interaction=nonstopmode -halt-on-error ms.tex > /dev/null; }

# 1. Generate ms.bbl from references.bib in a scratch copy.
B=$(mktemp -d); cp -r "$OUT"/. "$B"; cp references.bib "$B"
(cd "$B" && latex3 && bibtex ms > /dev/null && latex3 && latex3)
cp "$B/ms.bbl" "$OUT/ms.bbl"

# 2. Clean-room test: only the bundle's files, compiled three times.
T=$(mktemp -d); cp -r "$OUT"/. "$T"
(cd "$T" && latex3 && latex3 && latex3)
undef=$(grep -c -E 'undefined|Citation .* undefined' "$T/ms.log" || true)
pages_bundle=$(grep -o 'Output written on ms.pdf ([0-9]* pages' "$T/ms.log" | grep -o '[0-9]* pages')
abs=$(grep -c -E '(^|[^A-Za-z\\])[A-Za-z]:/|/Users/|/home/' "$T/ms.tex" || true)   # drive paths, not set-builder colons
tar -czf arxiv_bundle.tar.gz -C "$OUT" .
echo "bundle files: $(cd "$OUT" && find . -type f | sort | tr '\n' ' ')"
echo "bundle size: $(du -k arxiv_bundle.tar.gz | cut -f1) KB; clean build: $pages_bundle; undefined: $undef; absolute paths: $abs"
rm -rf "$B" "$T"
[ "$undef" -eq 0 ] && [ "$abs" -eq 0 ]
