#!/bin/sh
# Renders each diagram to docs/images as an SVG with glyphs as paths.
set -eu
cd "$(dirname "$0")"
out=../../../docs/images
tmp=$(mktemp -d)
for tex in *.tex; do
  latex -interaction=nonstopmode -halt-on-error -output-directory="$tmp" "$tex" >/dev/null
  dvisvgm --no-fonts --exact-bbox --output="$out/${tex%.tex}.svg" "$tmp/${tex%.tex}.dvi" 2>/dev/null
done
rm -rf "$tmp"
